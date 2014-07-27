Illumina_tools
==============

Tools for handling Illumina 

##Scripts/Programs Included

* DesignPrimers.py

###DesignPrimers.py

Helps with gene specific portion of the primer deisgn for Illumina paired end projects.

###Requirements:

 * Biopython (1.44 or newer - http://biopython.org/wiki/Download)

 * eprimer32 from the EMBOSS suite (http://emboss.sourceforge.net/download/)

 * primer32_core and primer3_config(http://sourceforge.net/projects/primer3/files/primer3/2.3.6/)

**Note- primer32_core is downloaded as primer3_core. You must change the name to primer32_core. This is to do to some wrapper functionallity from eprimer32.** 

**It's easiest to just put eprimer32, primer32_core, and primer3_config in the same directory that you call DesignPrimers.py in...I could write some copying and temp file functionality, but that is a lot of extra-uneeded work.**

###Examples

```
->./DesignPrimers.py --help
usage: Illumina Primer Design [-h] -fr full_region.fasta -ri
                              region_of_interest.fasta [-pr {150,300,500}]
                              [-in INSERTION]

A Program to Design Illumina Primers on the HiSeq, MiSeq, and NextSeq
Platforms

optional arguments:
  -h, --help            show this help message and exit
  -fr full_region.fasta, --full_region full_region.fasta
                        The full gene that we will PCR, ex. A whole antibody
                        region
  -ri region_of_interest.fasta, --region_of_interest region_of_interest.fasta
                        The region of interest that needs to be sequenced, ex.
                        The HCDR3 region
  -pr {150,300,500}, --paired_read_lengths {150,300,500}
                        The read length of the run. We automatically do the
                        paired end math for you
  -in INSERTION, --insertion INSERTION
                        Will there be any unexpected inseritons. For instance,
                        will the readlength of interest be any longer than
                        what you have specified?
```

1. You have the full sequence and the region of interest. The full sequence would be a vector, antibody, or extended portion of the genes you will design the primers on. The sequence of interest is the portin you wish to be resolved in the high throughput sequencing. Each should be contained in its own FASTA file.

```

>Full_gene
CAGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCAATG
GAACCAGCAATGATGTTGGTGGCTATGAATCTGTCTCCTGGTACCAACAACATCCCGGCAAAGCCCCCAA
AGTCGTGATTTATGATGTCAGTAAACGGCCCTCAGGGGTTTCTAATCGCTTCTCTGGCTCCAAGTCCGGC
AACACGGCCTCCCTGACCATCTCTGGGCTCCAGGCTGAGGACGAGGGTGACTATTACTGCAAGTCTCTGA
CAAGCACGAGACGTCGGGTTTTCGGCACTGGGACCAAGCTGACCGTTCTA

>Portion_you_care_about
ATGAATCTGTCTCCTGGTACCAACAACATCCCGGCAAAGCCCCCAAAGTCGTGATTTATGATGTCAGTAAA
CGGCCCTCAGGGGTTTCTAATCGCTTCTCTGGCTCCAAGTCCGGCAACACGGCCTCCCTGACCATCTCTGG
```

Then you can run the command:

```
python DesignPrimers.py -ri test/portion_you_care_about.fasta -fr test/full_gene.fasta

Your run parameters:
--------------------
Full Region:
AGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCAATGGAACCAGCAATGATGTTGGTGGCTATGAATCTGTCTCCTGGTACCAACAACATCCCGGCAAAGCCCCCAAAGTCGTGATTTATGATGTCAGTAAACGGCCCTCAGGGGTTTCTAATCGCTTCTCTGGCTCCAAGTCCGGCAACACGGCCTCCCTGACCATCTCTGGGCTCCAGGCTGAGGACGAGGGTGACTATTACTGCAAGTCTCTGACAAGCACGAGACGTCGGGTTTTCGGCACTGGGACCAAGCTGACCGTTCTA

Region of Interest:
ATGAATCTGTCTCCTGGTACCAACAACATCCCGGCAAAGCCCCCAAAGTCGTGATTTATGATGTCAGTAAACGGCCCTCAGGGGTTTCTAATCGCTTCTCTGGCTCCAAGTCCGGCAACACGGCCTCCCTGACCATCTCTGG

Posible Insertion Lenght:0
Read Length:300
Paird-End Read Length:150
-----------------------

Rank: 1
---------
Forward Primer: TCTGGGTCTCCTGGACAGTC
Forward Primer Length: 20
Forward GC Content: 60.0
Forward Primer Tm: 60.25
Reverse Primer: GACGTCTCGTGCTTGTCAGA
Reverse Primer Length: 20
Reverse GC Content: 55.0
Reverse Primer Tm: 60.04
Forward Region of Interest Coverage: 86
Reverse Region of Interest Coverage: 91
Minimum Paired-End Overlap: 35
Maximum Paired-End Overlap: 35

```

2. If you had a smaller region of interest you could use the 150 kit on the NextSeq with the -pr

```
>Full_gene
CAGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCAATG
GAACCAGCAATGATGTTGGTGGCTATGAATCTGTCTCCTGGTACCAACAACATCCCGGCAAAGCCCCCAA
AGTCGTGATTTATGATGTCAGTAAACGGCCCTCAGGGGTTTCTAATCGCTTCTCTGGCTCCAAGTCCGGC
AACACGGCCTCCCTGACCATCTCTGGGCTCCAGGCTGAGGACGAGGGTGACTATTACTGCAAGTCTCTGA
CAAGCACGAGACGTCGGGTTTTCGGCACTGGGACCAAGCTGACCGTTCTA

>Portion_you_care_about_smaller
ATGAATCTGTCTCCTGGTACCAACAACATCCCGGCAAAGCCCCCAAAGTCGTGATTTATGATGTCAGTAAA
CGG

>
python DesignPrimers.py -ri test/portion_you_care_about_smaller.fasta -fr test/full_gene.fasta -pr 150

Rank: 1
---------
Forward Primer: TCTGGGTCTCCTGGACAGTC
Forward Primer Length: 20
Forward GC Content: 60.0
Forward Primer Tm: 60.25
Reverse Primer: AATCACGACTTTGGGGGCTT
Reverse Primer Length: 20
Reverse GC Content: 50.0
Reverse Primer Tm: 59.89
Forward Region of Interest Coverage: 11
Reverse Region of Interest Coverage: 93
Minimum Paired-End Overlap: 30
Maximum Paired-End Overlap: 30

```

3. However, say your area of interest can have insertions. You can mark that with -in flag.


```

python DesignPrimers.py -ri test/portion_you_care_about_smaller.fasta -fr test/full_gene.fasta -in 30

Rank: 1
---------
Forward Primer: TGGTACCAACAACATCCCGG
Forward Primer Length: 20
Forward GC Content: 55.0
Forward Primer Tm: 59.96
Reverse Primer: CGGACTTGGAGCCAGAGAAG
Reverse Primer Length: 20
Reverse GC Content: 60.0
Reverse Primer Tm: 60.11
Forward Region of Interest Coverage: 89
Reverse Region of Interest Coverage: 35
Minimum Paired-End Overlap: 50
Maximum Paired-End Overlap: 80


```

What do these mean.

Most are self explanatory, but the coverage is how many bases into the foward region of interest will get. The same for the reverse.

The overlap is how much the paired in reads will overlap each other in the region of interest. There is a minimum and maximum in case there are insertions. The insertions will make the minimum less, but the program will complensate. If you don't like the overlap, try changing the read length kit.


