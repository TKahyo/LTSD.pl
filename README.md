#LTSD.pl
  
##Description  
#####LTSD.pl is a Perl script to estimate putative target site duplications (TSDs) of LTR-retrotransposons. 
  
##Requirement
#####1. Tools  
* BEDtools (flankBed and fastaFromBed): [Quinlan Lab at the University of Utah](http://bedtools.readthedocs.io/en/latest/)	v2.25.0 or later.  
* T-Coffee: [Notredame Lab, Comparative Bioinformatics Group at the Bioinformatics and Genomics Programme Center for Genomic Regulation](http://www.tcoffee.org/Projects/tcoffee/#Download) Version 11.00.8cbe486 or later  

#####2. Data sets  
* Reference genome (genome.fa): Downloadable from [the University of California Santa Cruz (UCSC) Genome Browser](http://genome.ucsc.edu/index.html)  
* Chromosome size data (genome\_num.txt): [chromName] [TAB] [chromSize]    

#  
    chr1    249250621  
    chr2    243199373  
        ...  

* ORF data (ORF\_list.bed): Annotation data of open reading frames (ORFs). Downloadable via RepeatMasker track from the Table Browser in  [the University of California Santa Cruz (UCSC) Genome Browser](http://genome.ucsc.edu/index.html). If LTR5_Hs elements are analyzed, HERVK-int should be selected.  

#  
    chr1    12840257    12845090    HERVK-int    27617    -  
    chr1    13459295    13460029    HERVK-int    4505     +  
        ...  

* LTR list data (LTR\_list.bed): [chromName] [TAB] [Start] [TAB] [End] [TAB] [+/-] [TAB] [insert]. If the LTRs absent on the reference genome are added in this list in order to investigate neighbor genes using Gene list data (see below), 'insert' should be replaced to 'pre' in the last column, and TSD genome positions should be written as following:  

#  
    chr1    1345186     1346153     +    insert  
    chr1    79792628    79792633    -    pre  
        ...  

* Gene list data (Gene\_list.bed, optional): [chromName] [TAB] [Start] [TAB] [End] [TAB] [GeneName] [TAB] [Strand].  

#  
    chr1    66999824    67210768    SGIP1    +  
    chr1    48998526    50489626    AGBL4    -  
        ...  

  
##Usage
    perl LTSD.pl -gf /path/genome.fa -gn /path/genome_num.txt -i /path/ORF_list.bed -tc /path/t_coffee -flank /path/flankBed -fasta /path/fastaFromBed [other_options] /path/LTR_list.bed  

#####Options
-h|--help　　print help  

-b　　Maximum of flanking sequence lengths (default 30)

*-flank|-fasta　　flankBed and fastaFromBed paths

*-gene　　File path of Gene\_list.bed

*-gf　　File path of genome.fasta

*-gn　　File path of genome\_num.txt

*-i　　File path of ORF\_list.bed

-s|-m　　Column numbers of strand and insert/pre in the LTR\_list.bed (default -s 4 -m 5)

-t　　Thread number (default -t 1)

*-tc　　T-Cofee path

*required options

##Demo
Demo datasets for LTR5_Hs on the human chromosome 19 (GRCh37/hg19) are included in the 'DEMO' directory.  

    perl LTSD.pl -b 100 -t 4 -gf chr19.fa -gn genome_num.txt -i ORF_list.bed -tc /path/t_coffee -flank /path/flankBed -fasta /path/fastaFromBed LTR5_Hs_chr19.bed

It took about 13 minutes under the condition of Intel i7-3930K CPU @ 3.20GHz multi cores and 64GB  memory.

Two result files are obtained.  

* tsd\_results.txt  
    Putative TSD sequences estimated from LTR\_list.bed.  
* tsd\_results_fil.txt  
    Putative TSD sequences filtered by a criteria.  
    
These files are included 'DEMO/results' directory.  

1st column: LTR\_inf | SOLO/PRE/paired-LTRs\_inf | TSD-insert-TSD\_inf | Present/Tsd | Identity.  
The other columns: [TSD\_seq]    [Left\_flanking\_seq]    [Right\_flanking\_seq]    [Distance\_to\_Gene]    [Gene]    [Direction]    [Distance\_to\_Gene]    [Gene]    [Direction]    [N\_count\_of\_TS].



