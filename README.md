Illumina_tools
==============

Tools for handling Illumina 

#Scripts Programs Included
==========================

* DesignPrimers.py

##DesignPrimers.py

Helps with gene specific portion of the primer deisgn for Illumina paired end projects.

#Requirements:

*Biopython (1.44 or newer - http://biopython.org/wiki/Download)
*eprimer32 from the EMBOSS suite (http://emboss.sourceforge.net/download/)
*primer32_core and primer3_config(http://sourceforge.net/projects/primer3/files/primer3/2.3.6/)

**Note primer32_core is downloaded as primer3_core. You must change the name to primer32_core. This is to do to some wrapper functionallity from eprimer32. 

**It's easiest to just put eprimer32, primer32_core, and primer3_config in the same directory that you call DesignPrimers.py in...I could write some copying and temp file functionality, but that is a lot of extra-uneeded work.

#Examples

'''
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
'''