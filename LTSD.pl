#!/usr/local/bin/perl

#####################################################
# Estimating putative TSDs from LTR-retrotransposon list.
# perl LTSD.pl -gf genome.fa -gn genome_num.txt -i HERVint.bed -tc t_coffee -flank flankBed -fasta fastaFromBed [other_options] LTR_list.bed
# Below tools and data are required in this script.
# [1] Reference genome (.fa)
# [2] bedtools
# [3] Genome file utilized for flankbed (<chromName><TAB><chromSize>)
# [4] T-Coffee
# [5] LTR list (.bed)	# see README
# [6] ORF list (.bed)	# from RepeatMasker
# [7] Gene list (.bed)	# optional, see README
#####################################################

use strict;
use warnings;
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);
use Pod::Usage;
use IO::Handle;
use threads;

STDOUT->autoflush(1);
my $start_time = (time())[0];
print $0, "\n";

# Default path setting
my $genome_fa;
my $flankbed;
my $fastaFromBed;
my $genome_num;
my $t_coffee;
my $int;
my $gene_name; # optional

my $th = 1;		# number of threads
my $s_col = 4; 		# strand column
my $m_col = 5;		# pre/insert column
my $flank_bp = 30; 	# used for flankbed
my %opt = ();
GetOptions(\%opt, 't=i' => \$th,'s=i' => \$s_col, 'm=i' => \$m_col, 'i=s' => \$int, 
'b=i' => \$flank_bp, 'gn=s' => \$genome_num, 'gf=s' => \$genome_fa, 'gene=s' => \$gene_name, 'tc=s' => \$t_coffee, 'flank=s' => \$flankbed, 'fasta=s' => \$fastaFromBed, 'help|h');
pod2usage(1) if ($opt{help} or $opt{h});

# OTHER INTERNAL PARAMETERS
my $too_long_distance = 12000;	# distance between two LTRs
my $thresh_LTR_length = 4;
my $gap_LTR_int_min = -1;
my $gap_LTR_int_max = 100;
my $thresh_score_proviral_TSD = 0.5;
my $iscore = 0.75;


my $bed_file = shift @ARGV;
die "ERROR: bed_files of LTR positions is required\n" if(!defined $bed_file);
die "ERROR: Optioins -gf (genome.fa) and -gn (genome position file) are required\n" if(!defined $genome_fa or !defined $genome_num);
die "ERROR: Option -i (e.g. HERVK-int.bed) is necessary\n" if(!defined $int);
print "THREADS: $th\n";
$s_col--;
$m_col--;


print "> Sorting LTR position data...\n";
my $bed_file_sorted = 'input_sorted.tmp';
system("sort -V -k 1,1 -k 2,2n $bed_file >$bed_file_sorted");


print "> Importing bed_data...\n";
my @data_th; # (thread_0_hash, thread_1_hash, ... ); hash->{$key} = \@content; 
my @list_th; # (thread_0_array, thread_1_array, ... ); array = (key1, key2, ...)
# e.g. LTR data (13 rows) => thread_0: 1~5, thread_1: 5=8, thread_3: 9~12 (in the case of $th=3)  
my $wc = `wc -l $bed_file_sorted | cut -f1 -d' '`;
chomp($wc);
my @nnn; # number list of distributed data
my $i = 0; # thread number
foreach(my $r=0; $r < $th; $r++) {
	$nnn[$r] = int($wc/$th);
}
my $rem = $wc - int($wc/$th)*$th;
while(1) {
	last if($rem == 0);
	$nnn[$i]++;
	$i++;
	$rem--;
}

my $row = 1;
$i = 0;
my @previous_pos = ('chr1',0,0,'+'); # initial value
open my $bb, '<', $bed_file_sorted or die $!;
while(<$bb>) {
	chomp;
	my @content = split /\t/, $_;
	my $key = join ":", @content;
	
	# redistribution of LTR data
	my $sum = 0;
	map {$sum+=$_} @nnn[0..$i];
	# considering proviral state
	if($row > $sum) {
		if($content[0] ne $previous_pos[0] or $content[$s_col] ne $previous_pos[3] or $content[1]-$previous_pos[2] > $too_long_distance) {
			my $pi = $i;
			$i++;
			my $pre_row = $row - 1;
			print "THE_LAST_LINE_OF_THREAD-$pi: $pre_row\n";
			print "BORDER_DATA [[$pi-$i]] >>\n";
			print join "\t", @previous_pos, "\n";
			print join "\t", @content[0..3],"\n";
		} elsif($nnn[$i+1] <= 1) {
			warn "REMAIN <= 1 \@ $th\n";
			$i++;
		} else {
			$nnn[$i]++;
			$nnn[$i+1]--;
		}
	}
	$data_th[$i]->{$key} = $_;
	push @{$list_th[$i]}, $key;
	
	@previous_pos = @content[0..2,$s_col];
	$row++;
}
close $bb;

print "Thread: number of distributed data\n";
foreach(my $n=0; $n < @nnn; $n++) {
	print "$n: $nnn[$n]\n"; 
}

open my $ii, '<', $int or die $!;
my $int_data; # hash of HERVK-int
while(<$ii>) {
	chomp;
	my @content = split /\t/, $_;
	my $key = join ":", @content;
	$int_data->{$key} = $_;
}
close $ii;


print "> flankBed & fastaFromBed...\n";
my $tmp1 = 'flank_out.tmp';
my $tmp2 = 'fast_out.tmp'; # flanking sequence
system("$flankbed -i $bed_file_sorted -g $genome_num -b $flank_bp >$tmp1");
print "=====\n";
system("$fastaFromBed -fi $genome_fa -bed $tmp1 -fo $tmp2");

my @tsd_data_th; # (thread_0_hash, thread_1_hash, ... ); hash->{LTR_pos} = flanking_seq[(LEFT,RIGHT)]
$i = 0;
my $line = 1;
my $thresh_index = 0; # array index of LTR list
my $seq;
my $sum = $nnn[$i];
open my $t2, '<', $tmp2 or die $!;
while(<$t2>) {
	if($i < $th and $line/4 > $sum) {
		my $pi = $i;
		$i++;
		print "$pi => $i\n";
		$sum += $nnn[$i];
		$thresh_index = 0; # reset
	}
	chomp;
	
	if($line%2 == 0) {
		$seq = $_;
	} else {
		$line++;
		next;
	}
	
	if(exists $list_th[$i]->[$thresh_index]){ # LTR_pos
		push @{$tsd_data_th[$i]->{$list_th[$i]->[$thresh_index]}}, $seq;
		$thresh_index++ if($line%4 == 0);
	} else {
		print $line/4, "\n";
		die "ERROR\@$i\n$line\n$thresh_index\n$sum\n";
	}
	$line++;
}
close $t2;


print "> TSD assesment & saving...\n";
my $save = 'tsd_results.txt';
my @threads;
for (my $i = 0; $i < $th; $i++) {
	$threads[$i] = threads->new(\&body, $list_th[$i], $tsd_data_th[$i], $i);

}

for (my $i = 0; $i < $th; $i ++) {
	my $rtn = $threads[$i]->join;
	warn "::THREAD_JOIN_ERROR\n" if ($rtn != 0);
}

system("cat tsd_out*.txt >$save");
system("rm tsd_out_*.txt");
system("rm *.tmp");


print "> Filtering ...\n";
(my $save1 = $save) =~ s/\.(.+)$/_N-count\.$1/;
(my $save2 = $save) =~ s/\.(.+)$/_fil\.$1/;

open my $out2, '>', $save2 or die $!;
open my $out, '>', $save1 or die $!;
open my $in, '<', $save or die $!;
while(<$in>){
	chomp;
	my $count = 0;
	my ($first_col, $tsd) = (split /\t/, $_)[0,1];
	my $score = (split /\|/, $first_col)[4]; 
	if(!defined $tsd) {
		print $_,"\n";
	}
	while($tsd=~/N/g){
		$count++;
	}
	print $out "$_\t$count\n";
	if($score >= $iscore and ((length($tsd) >=4 and length($tsd) <= 10 and $count <=1) or (length($tsd) >=11 and length($tsd) <=20 and $count <=2) or (length($tsd) >=21 and $count <=3)) and $tsd !~ /provirus/) {
		print $out2 "$_\t$count\n";
	}
}
close $in;
close $out;
close $out2;


system("mv $save1 $save");
print "> 1st column: LTR_inf|SOLO/PRE/paired-LTRs_inf|TSD-insert-TSD_inf|Present/Tsd|Score\n";
print "> The other columns: TSD_seq\tLeft_flanking_seq\tRight_flanking_seq\tDistance_to_Gene\tGene\tDirection\tDistance_to_Gene\tGene\tDirection\tN_count_of_TSD\n";
print "> OUTPUT: tsd_results.txt and tsd_results_fil.txt\n";

#======================================================================

sub body {
	my ($l, $tsd_data, $i) = @_;
	my $solo_flag = 1;
	my @list = @{$l}; # LTR_pos list
	my $save_body = 'tsd_out_'.$i.'.txt';
	my $output_tsd = 'output_tsd_'.$i.'.tmp';
	system("touch $output_tsd");
	my $input_tsd = 'input_tsd_'.$i.'.tmp';
	my $ncopy;
	
	my $dir_i = 'aln_'.$i;
	system("mkdir $dir_i");
	
	open my $o, '>', $save_body or die $!;
	foreach (my $n = 0; $n < @list; $n++) {
		$ncopy = $n;
		print $o "$list[$n]";
		
		if($solo_flag <= 0) {
			my $allele = "provirus";
			$allele = "tandem" if($solo_flag < 0);
			print $o "\t$allele\n";
			$solo_flag++;
			next;
		}
		
		my ($chr, $start, $end, $strand, $mark) = (split /:/, $list[$n])[0..2, $s_col, $m_col];
		my $tsd_seqL = $tsd_data->{$list[$n]}[0];
		my $tsd_seqR = $tsd_data->{$list[$n]}[1];
		my $tsd; # common sequence
		my ($tsd_1, $tsd_2) = ("ND", "ND"); # Left and Right flanking sequences
		
		# LTR status (Present/Short/Tsd)
		my $hg19 = 'Present';
		my $len = $end - $start;
		$hg19 = 'Short' if($len < $thresh_LTR_length);
		$hg19 = 'Tsd' if($mark =~ /pre/);
		
		# Addressing HERVint
		# 0: ND, 1:LTR(+)-HERVint(+), 2: HERVint(+)-LTR(+), 3: HERVint(-)-LTR(+), 4: HERVint(-)-LTR(-)
		my $ltr_side = 0;
		for my $i (keys %{$int_data}) {
			my ($c, $s, $e, $d) = (split /:/, $i)[0..2,5];
			next if($chr ne $c or $strand ne $d);
			if($d eq '+') {
				$ltr_side = 1, last if($s - $end >= $gap_LTR_int_min and $s - $end <= $gap_LTR_int_max);
				$ltr_side = 2, last if($start - $e >= $gap_LTR_int_min and $start - $e <= $gap_LTR_int_max);
			} elsif($d eq '-') {
				$ltr_side = 3, last if($s - $end >= $gap_LTR_int_min and $s - $end <= $gap_LTR_int_max);
				$ltr_side = 4, last if($start - $e >= $gap_LTR_int_min and $start - $e <= $gap_LTR_int_max);
			}
		}
		
		# defining flanking sequences
		my $add = 1; # solo: 0, putative provirus: 1, putative tandem: >=2
		my $first = $tsd_seqL;
		my $second = $tsd_seqR;
		my ($chr2, $start2, $end2, $strand2);
		my $end3;
		my $score;
		
		open my $it, '>', $input_tsd or die $!;
		while(1) {
			if($hg19 eq 'Short') {
				$tsd = 'UC'; # Unconfirmed
				$add--;
				last;
			} elsif($hg19 eq 'Tsd') {
				my $tsd_out = join "\t", ($chr, $start, $end, $strand);
				print $it $tsd_out, "\n";
				$add--;
				last;
			}
			
			# SOLO
			if(!exists $list[$n+$add] or $strand !~ /[-+]/ or $ltr_side == 2 or $ltr_side == 4) {
				($tsd, $score, $tsd_1, $tsd_2) = @{&tsd($first, $second, \@list, $ncopy, $dir_i)};
				#print "A: Side=$ltr_side => $tsd_1\t$tsd_2\n";
				$add--;
				last;
			}	
			
			# partner_LTR
			($chr2, $start2, $end2, $strand2) = (split /[:]/, $list[$n+$add])[0..2,$s_col];
			my $distance_partner = $start2 - $end;
			my $partner_index = 1;
			if($strand2 !~ /[-+]/) {
				($tsd, $score, $tsd_1, $tsd_2) = @{&tsd($first, $second, \@list, $ncopy, $dir_i)};
				#print "A: Side=$ltr_side => $tsd_1\t$tsd_2\n";
				$add--;
				last;
			}
			my $tsd_next = $tsd_data->{$list[$n+$add]}[$partner_index]; # right flanking sequence of partner_LTR
			my $decision;
			($decision, $score, $tsd_1, $tsd_2) = @{&tsd($first, $second, \@list, $ncopy, $dir_i)};
			# if($tsd_1 =~ /ND/) {
				# print ":::ND\n";
			# }
			$add--, $tsd = $decision, last if($chr ne $chr2 or $strand ne $strand2 or $distance_partner > $too_long_distance*$add); # SOLO
			my ($decision_next, $score_next, $tsd_1_next, $tsd_2_next) = @{&tsd($first, $tsd_next, \@list, $ncopy, $dir_i)}; # Probably provirus	
			
			## final decision for SOLO or provirus
			if($score > $score_next or $score_next < $thresh_score_proviral_TSD) { # SOLO
				$add--;
				$tsd = $decision;
				last;
			}
			$solo_flag--;
			$second = $tsd_next;
			$tsd = $decision_next;
			$tsd_1 = $tsd_1_next;
			$tsd_2 = $tsd_2_next;
			$add++;
			$end3 = $end2;
		}
		close $it;
		$tsd_1 = uc($tsd_1) if(defined $tsd_1);
		$tsd_2 = uc($tsd_2) if(defined $tsd_2);
		
		# pre on the reference genome
		if(!defined $tsd) {
			system("$fastaFromBed -fi $genome_fa -bed $input_tsd -fo $output_tsd");
			open my $ot, '<', $output_tsd or die $!;
			while(<$ot>) {
				chomp;
				next if($_ =~ /^>/);
				$tsd = uc($_);
			}
			close $ot;
			print "::TSD described in the bed file >> $tsd @ $list[$n]\n";
		}
		
		# TSD_position
		my ($tsd_chr, $tsd_s, $tsd_e) = ('-', '-', '-');
		if($solo_flag == 1) {
			my $status = 'SOLO';
			if($tsd ne 'ND' and $tsd ne 'UC') {
				$tsd_chr = $chr;
				$tsd_s = $start - length($tsd);
				$tsd_e = $end + length($tsd);
			}
			if($hg19 eq 'Tsd') {
				$tsd_s = $start;
				$tsd_e = $end;
				$status = 'PRE';
			}
			my $tr = join ":", ($tsd_chr, $tsd_s, $tsd_e, $strand);
			if($tsd ne 'ND' and $strand eq '-') {
				$tsd =~ tr/[ATGC]/[TACG]/;
				$tsd = reverse($tsd);
				if($tsd_1 ne 'ND') {
					$tsd_1 =~ tr/[ATGC]/[TACG]/;
					$tsd_1 = reverse($tsd_1);
				}
				if($tsd_2 ne 'ND') {
					$tsd_2 =~ tr/[ATGC]/[TACG]/;
					$tsd_2 = reverse($tsd_2); 
				}
			}
			
			## Neighbor gene
			my $distance_inf = "-\t-\t-\t-";
			if(defined $gene_name) {
				$distance_inf = join "\t", @{&distance($gene_name, $tr, $ltr_side)};
			}
			print $o "|$status|$tr|$hg19|$score\t$tsd\t$tsd_1\t$tsd_2\t$distance_inf\n";
		} elsif($solo_flag != 1) {
			if($tsd ne 'ND' and $tsd ne 'UC') {
				$tsd_chr = $chr;
				$tsd_s = $start - length($tsd);
				$tsd_e = $end3 + length($tsd);
			}
			my $tr = join ":", ($tsd_chr, $tsd_s, $tsd_e, $strand);
			if($tsd ne 'ND' and $strand eq '-') {
				$tsd =~ tr/[ATGC]/[TACG]/;
				$tsd = reverse($tsd);
				if($tsd_1 ne 'ND') {
					$tsd_1 =~ tr/[ATGC]/[TACG]/;
					$tsd_1 = reverse($tsd_1);
				}
				if($tsd_2 ne 'ND') {
					$tsd_2 =~ tr/[ATGC]/[TACG]/;
					$tsd_2 = reverse($tsd_2); 
				}
			}
			
			## Neighbor gene
			my $distance_inf = "-\t-\t-\t-";
			if(defined $gene_name) {
				$distance_inf = join "\t", @{&distance($gene_name, $tr, $ltr_side)};
			}
			print $o "|$list[$n+$add]|$tr|$hg19|$score\t$tsd\t$tsd_1\t$tsd_2\t$distance_inf\n";
		}
	}
	close $o;
	system("rm -r $dir_i");
}

sub tsd {
	my $t1 = shift; # LEFT_LTR
	my $t2 = shift; # RIGHT_LTR
	my $list = shift;
	my $ncopy = shift;
	my $dir_i = shift;
	my $hash; # key:$n, value:$max_i
	my $hash_max_seq; # key:$n, value:[($max_seq1, $max_seq2)]
	
	die "ERROR: Different length bwtween $t1 (ts1) and $t2 (ts2)\n" if(length($t1) != length($t2));
	my $max = -20; # minimum
	#print "$t1 <LTR> $t2\n";
	
	# calculating alignment score (tsd_length = n)
	foreach(my $n = length($t1); $n > 0; $n--) {
		## Alignment between t1(-2,-1,0,+1,+2) and t2 using T-COFFEE
		## [$n:4~7] -1,+1
		## [$n:8~(len-3)] -2,-1,+1,+2
		## [$n:len-2] +1,-1,-2
		## [$n:len-1] -1,-2
		my @t1_right_nA = (substr($t1, (-1)*$n));
		if($n >=4 and $n <=7) {
			unshift @t1_right_nA, substr($t1, (-1)*$n+1);
			push @t1_right_nA, substr($t1, (-1)*$n-1);
		} elsif($n >=8 and $n <= length($t1)-3) {
			unshift @t1_right_nA, (substr($t1, (-1)*$n+2), substr($t1, (-1)*$n+1));
			push @t1_right_nA, (substr($t1, (-1)*$n-1), substr($t1, (-1)*$n-2));
		} elsif($n == length($t1)-2) {
			unshift @t1_right_nA, (substr($t1, (-1)*$n+2), substr($t1, (-1)*$n+1));
			push @t1_right_nA, (substr($t1, (-1)*$n-1));
		} elsif($n == length($t1)-1) {
			unshift @t1_right_nA, (substr($t1, (-1)*$n+2), substr($t1, (-1)*$n+1)); #
		}
		my $t2_left_n = substr($t2, 0, $n+1);
		
		my $alined;
		my @get_score;
		my $max_i = -50*length($t1); # minimum	
		foreach(my $i=0; $i < @t1_right_nA; $i++) {
			my ($name_1, $name_2) = ("t1_".$n.'_'.$i, "t2_".$n.'_'.$i);
			$alined = &tcoffee($name_1, $t1_right_nA[$i], $name_2, $t2_left_n, $alined, $dir_i);
			$get_score[$i] = &scoring($alined->{$name_1}, $alined->{$name_2});
			#print ":::$i=>$get_score[$i]\n[1]$alined->{$name_1}\n[2]$alined->{$name_2}\n\n";
			$max_i = $get_score[$i] if($max_i < $get_score[$i]);
		}		
		#print "$n>>$max_i\n1:$t1_right_nA[$max]\n2:$t2_left_n\n\n";
		
		## the longest sequence of $max_i
		my $max_seq1 = "ND";
		my $max_seq2 = "ND";
		foreach(my $ii = 0; $ii <=$#get_score; $ii++) {
			if($get_score[$ii] == $max_i) {
				my $key1 = "t1_".$n.'_'.$ii;
				my $key2 = "t2_".$n.'_'.$ii;
				
				if(exists $alined->{$key1} and exists $alined->{$key2}) {
					$max_seq1 = $alined->{$key1} ;
					$max_seq2 = $alined->{$key2};	
				} else {
					warn "*** WRANING: $key1/$key2 values not exist ***\n";
				}
			}
		}
		
		if($max_seq1 eq 'ND') {
			my @a = values %{$alined};
			warn "==== ND ==== $max_i <=> @get_score\t@a \n";
			### if this warning is shown, coordinate the value of $max_i
		}
		
		push @{$hash_max_seq->{$n}}, ($max_seq1, $max_seq2);
		## processing in the case of $n<=3 && $max_1==1
		if($n <= 3 and $max_i == 1 and $max >= $iscore) {
			print ">> SELECT LONG TSD: score $max => $list->[$ncopy]\n";
			last;
		}
		$max = $max_i if($max < $max_i);
		$hash->{$n} = $max_i;
		
		last if($max_i == 1);
		
		$max_i = 0;
	}
	# for my $k (sort {$a <=> $b} keys %{$hash}) {
		# print "$k\t$hash->{$k}\t$max\n";
	# }
	
	my @max_grep = grep {$hash->{$_} == $max} keys %{$hash};
	my @sorted = sort {$b <=> $a} @max_grep; # longest TSD
	
	my $r1 = $hash_max_seq->{$sorted[0]}[0];
	my $r2 = $hash_max_seq->{$sorted[0]}[1];
	my $r3 = uc($r1);
	$r3 = uc($r2) if(length($r1) < length($r2));
	## replace variants to 'N'
	if($r1 ne 'ND') {
		if(defined $r1 and $r1 !~ /^$r2$/i) {
			foreach(my $l = 0; $l < length($r3); $l++) {
				my $a = substr($r1, $l, 1);
				my $b = substr($r2, $l, 1);
				if($a !~ /^$b$/i) {
					substr($r3, $l, 1) = "N";
				}
			}
		}
	}
	
	return [($r3, $max, $r1, $r2)];
}

sub scoring {
	my ($t1, $t2) = @_;
	my $score = 0;
	
	my $longer = length($t2);
	$longer = length($t1) if(length($t1) > length($t2));
	foreach(my $j = 0; $j < $longer; $j++) {
		my $m = - ($longer - $j);
		my $t1_base = substr($t1, $m, 1); 
		my $t2_base = substr($t2, $j, 1); 
		# print "$t1_base <=> $t2_base\n";
		# print "$t1\n$t2\n\n";
		if($t1_base =~ /^$t2_base$/i) { ;
			$score++;
		#} elsif($t1_base !~ /-/ and $t2_base !~ /-/) {
		#	$score += 0;
		} elsif($t1_base =~ /-/ or $t2_base =~ /-/) {
			$score -= 1; # gap penalty
		}
	}
	(my $raw_seq = $t2) =~ s/-//g;
	my $raw_len = length($raw_seq);
	$score = $score/$raw_len;
	return $score;
}

sub tcoffee {
	my ($seq1_name, $seq1, $seq2_name, $seq2, $inhash, $dir_i) = @_;
	
	# making .fa
	my $tmp_fa = $dir_i.'/'.$seq1_name.'_tmp_tcoffee.fa';
	open my $tc, '>', $tmp_fa or die $!;
	my $tname1 = '>'.$seq1_name;
	my $tname2 = '>'.$seq2_name;
	print $tc "$tname1\n$seq1\n$tname2\n$seq2\n";
	close $tc;
	
	# t-coffee
	my $tmp_log = $dir_i.'/'.$seq1_name.'_tmp_log';
	my $tmp_aln = $dir_i.'/'.$seq1_name.'_tmp_tcoffee.fasta_aln';
	system("$t_coffee $tmp_fa -type=dna -output fasta_aln -outfile $tmp_aln >$tmp_log 2>&1");
	while(1) {
		if(!-f $tmp_aln) {
			sleep 0.1;
		} else {
			last;
		}
	}
	unlink $tmp_log;
	my @get_alns = ();
	my $aln_flag = 0;
	open my $fa, '<', $tmp_aln or die $!;
	while(<$fa>){
		$aln_flag++, next if($_=~/^>/);
		chomp;
		$get_alns[$aln_flag] .= $_;
	}
	close $fa;
	
	# Output
	$inhash->{$seq1_name} = $get_alns[1];
	$inhash->{$seq2_name} = $get_alns[2];
	
	unlink $tmp_fa;
	unlink $tmp_aln;

	return $inhash;
}

sub distance {
	my $gene_inf = shift;
	my $region = shift;
	my $ltr_s = shift;
	
	my $ltr_up_distance;
	my $ltr_down_distance;
	my ($chr, $pos1, $pos2) = (split /:/, $region);
	$chr = $chr . '[[:space:]]' if($chr =~ /chr[\d\w]+$/);
	if($ltr_s eq '+') {
		chomp($ltr_up_distance = `grep $chr $gene_inf | perl -lane 'print abs($pos1-\$F[2]), "\t\$F[3]\t\$F[4]"' | sort -k1 -g | head -1`);
		chomp($ltr_down_distance = `grep $chr $gene_inf | perl -lane 'print abs(\$F[1]-$pos2), "\t\$F[3]\t\$F[4]"' | sort -k1 -g | head -1`);
	} else {
		chomp($ltr_down_distance = `grep $chr $gene_inf | perl -lane 'print abs($pos1-\$F[2]), "\t\$F[3]\t\$F[4]"' | sort -k1 -g | head -1`);
		chomp($ltr_up_distance = `grep $chr $gene_inf | perl -lane 'print abs(\$F[1]-$pos2), "\t\$F[3]\t\$F[4]"' | sort -k1 -g | head -1`);
	}

	return [($ltr_up_distance, $ltr_down_distance)];
}

print "done\n";

my $end_time = (time())[0];
my $process_time = $end_time - $start_time;
printf("Process time [%02d:%02d:%02d]\n", int($process_time/3600), int(($process_time%3600)/60), $process_time%60);

exit;


__END__
 
=pod
 
=head1 SYNOPSIS
 
B<perl LTSD.pl -gf /path/genome.fa -gn /path/genome_num.txt -i /path/ORF_list.bed -tc /path/t_coffee -flank /path/flankBed -fasta /path/fastaFromBed [other_options] /path/LTR_list.bed>

Options: [-h|--help|-b|-gene|-m|-p|-s|-t|-gf*|-gn*|i*|-tc*|-flank*|-fasta*; *required]

=head1 OPTIONS
 
=over 8
 
=item B<-h|--help>
 
print help.

=item B<-b>

Maximum of flanking sequence lengths.  (default 30)

=item B<-flank|-fasta>

flankBed and fastaFromBed paths

=item B<-gene>

File path of gene.bed. <chr> <start> <end> <geneName> <strand>

=item B<-gf>

File path of genome.fa.

=item B<-gn>

File path of genome.num, in which <chromosome> <start> <end> are described.

=item B<-i>

File path of ORF_list.bed

=item B<-s|-m>

Column numbers of strand and insert/pre in LTR_list.bed. (default -s 4 -m 5)

=item B<-t>

Thread. (default -t 1)

=item B<-tc>

T-COFFEE path

=back

=cut
