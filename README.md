#LTSD.pl
  
##Description  
#####LTSD.pl is a Perl script to estimate putative target site duplications (TSDs) of LTR-retrotransposons. 
  
##Requirement
### Tools  
 [1] BEDtools (flankBed and fastaFromBed): [Quinlan Lab at the University of Utah](http://bedtools.readthedocs.io/en/latest/)	v2.25.0 or later.  
 [2] T-Coffee: [Notredame Lab, Comparative Bioinformatics Group at the Bioinformatics and Genomics Programme Center for Genomic Regulation](http://www.tcoffee.org/Projects/tcoffee/#Download) Version 11.00.8cbe486 or later.  

### Data sets  
 [1] Reference genome (genome.fa): Downloadable from [the University of California Santa Cruz (UCSC) Genome Browser](http://genome.ucsc.edu/index.html)  
 [2] Chromosome size data (genome\_num.txt): [chromName] [TAB] [chromSize]    

#  
    chr1    249250621  
    chr2    243199373  
        ...  

 [3] Internal proviral sequence data (int\_list.bed): Annotation data of the internal provirals. Downloadable via RepeatMasker track from the Table Browser in  [the University of California Santa Cruz (UCSC) Genome Browser](http://genome.ucsc.edu/index.html). In the case of HML-2, HERVK-int should be included. The genomic positions are based on a zero-based start.  

#  
    chr1    12840257    12845090    HERVK-int    27617    -  
    chr1    13459295    13460029    HERVK-int    4505     +  
        ...  

 [4] LTR list data (LTR\_list.bed): [chromName] [TAB] [Start] [TAB] [End] [TAB] [+/-] [TAB] [insert]. In the case of HML-2, LTR5_Hs data should be included. The genomic positions are based on a zero-based start. If the LTRs absent on the reference genome are added in this list in order to investigate neighbor genes using Gene list data (see below), 'insert' should be replaced to 'pre' in the last column, and TSD genome positions should be written as following:  

#  
    chr1    1345186     1346153     +    insert  
    chr1    79792628    79792633    -    pre  
        ...  

##Demo
Demo datasets for LTR5_Hs/LTR5 on the human chromosome 19 (GRCh37/hg19) are included in the 'DEMO' directory except chr19.fa (58MB). LTR5_Hs_619.bed includes known LTR5_Hs/LTR5 data on autosomes.  

    perl LTSD.pl -b 100 -t 4 -gf chr19.fa -gn hg19_genome_num.txt -i int_list.bed -tc /path/t_coffee -flank /path/flankBed -fasta /path/fastaFromBed LTR5_Hs_chr19.bed

It took about 12 minutes under the condition of Intel i7-3930K CPU @ 3.20GHz multi cores and 64GB (DDR3 1333MHz) memory.

Two result files are generated.  

* tsd\_results.txt  
    Putative TSD sequences estimated from LTR\_list.bed.  
* tsd\_results_fil.txt  
    Putative TSD sequences filtered by a criteria.  
    
These files are included 'DEMO/results' directory.  
1st column: LTR\_inf | SOLO/PRE/paired-LTRs\_inf | TSD-insert-TSD\_inf | Present/Tsd | Score.  
The other columns: [TSD\_seq]    [Left\_flanking\_seq]    [Right\_flanking\_seq]    [Distance\_to\_Gene]    [Gene]    [Direction]    [Distance\_to\_Gene]    [Gene]    [Direction]    [N\_count\_of\_TSD].


##Usage
    perl LTSD.pl -gf /path/genome.fa -gn /path/genome_num.txt -i /path/int_list.bed -tc /path/t_coffee -flank /path/flankBed -fasta /path/fastaFromBed [other_options] /path/LTR_list.bed  

#####Options
-h|--help　　print help  

-b　　Maximum of flanking sequence lengths (default 30)  
　　　Minimum is fixed at 4.

-flank\*|-fasta\*　　/path/flankBed (-flank) and /path/fastaFromBed (-fasta) 

-gene　　File path of Gene\_list.bed to annotate neighbor genes  
　　　　　[chromName] [TAB] [Start] [TAB] [End] [TAB] [GeneName] [TAB] [Strand].

#  
    chr1    66999824    67210768    SGIP1    +  
    chr1    48998526    50489626    AGBL4    -  
        ...  

-gf\*　　genome.fasta

-gn\*　　genome\_num.txt

-i\*　　int\_list.bed

-s|-m　　Column numbers of strand and insert/pre in the LTR\_list.bed (default -s 4 -m 5)

-t　　Thread number (default -t 1)

-tc\*　　/path/T-Cofee path

*required
