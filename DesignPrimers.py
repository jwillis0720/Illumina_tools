import argparse
import tempfile
import sys
from collections import OrderedDict
try:
    import Bio.SeqIO
    from Bio.Emboss.Applications import Primer3Commandline as P3CL
    from Bio.Emboss import Primer3 as P3
    from Bio import pairwise2 as pw2
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
except ImportError:
    raise ImportError('Need BioPython to Run this Script')


class ArgumentParsingIllumina():

    """Simple Argument Parsing Class"""

    def __init__(self):
        self.arg_parse = argparse.ArgumentParser(
            prog="Illumina Primer Design",
            description="A Program to Design Illumina \
            Primers on the HiSeq, MiSeq, and NextSeq \
            Platforms")

        self.arg_parse.add_argument(
            '-fr', '--full_region', metavar="full_region.fasta",
            required=True,
            type=self.validate_single_fasta,
            help="The full gene that we will PCR, ex. A whole antibody region")

        self.arg_parse.add_argument(
            '-ri', '--region_of_interest', metavar='region_of_interest.fasta',
            type=self.validate_single_fasta,
            required=True,
            help="The region of interest that needs to be sequenced, ex. \
            The HCDR3 region")

        self.arg_parse.add_argument(
            '-pr', '--paired_read_lengths', choices=[150, 300, 500],
            type=int, default=300,
            help="The read length of the run. We automatically do the paired \
            end math for you")

        self.arg_parse.add_argument(
            '-in', '--insertion', type=int, default=0,
            help="Will there be any unexpected inseritons. For instance, \
            will the readlength of interest\
            be any longer than what you have specified?")

        # self.arg_parse.add_argument(
        #     '-st', '--stagger', default=False,
        #     action="store_true",
        #     help="Do you want to stagger primers to add ambiguous nucleotides \
        #     to the primers so they are not in the same frame")

        self.arg_parse = self.arg_parse.parse_args()
        self.validate_region_of_interest()

    def validate_single_fasta(self, file):
        try:
            Bio.SeqIO.read(open(file), 'fasta')
            return file
        except:
            raise argparse.ArgumentTypeError(
                '\n\nERROR:\nFile {0} is not a single fasta file'.format(file))

    def validate_region_of_interest(self):
        full_region = Bio.SeqIO.read(
            open(self.arg_parse.full_region), 'fasta').seq
        interest = Bio.SeqIO.read(
            open(self.arg_parse.region_of_interest), 'fasta').seq
        if interest not in full_region:
            raise argparse.ArgumentTypeError(
                "Region of Interest not in Full Region")


class IlluminaPrimer():

    def __init__(self):
        # command line arguments
        cla_dict = ArgumentParsingIllumina().arg_parse.__dict__

        # from argument dict
        self.paired_read_lengths = cla_dict['paired_read_lengths']
        self.full_region_file_name = cla_dict['full_region']
        self.region_of_interest_file_name = cla_dict['region_of_interest']
        self.insertion = cla_dict['insertion']
        self.region_length = cla_dict['paired_read_lengths']
        #self.stagger = cla_dict['stagger']

        # Needed for illumina
        # Primer needs to include this sequence
        self.seq_to_include = str(Bio.SeqIO.read(
            open(self.region_of_interest_file_name), 'fasta').seq).upper()
        self.full_seq = str(Bio.SeqIO.read(
            open(self.full_region_file_name), 'fasta').seq).upper()

        # run primer 3
        self.parse_primer3()

    def setup_primer3(self):
        setup_file = tempfile.NamedTemporaryFile()
        p3cl = P3CL(cmd='eprimer32')
        p3cl.sequence = self.full_region_file_name
        p3cl.outfile = setup_file.name
        p3cl.includedregion = self.seq_to_include
        p3cl.psizeopt = self.region_length - (self.insertion + 30)
        stdout, stderr = p3cl()
        return setup_file

    def parse_primer3(self):
        self.p3 = P3.read(self.setup_primer3()).primers

    def print_all_primers_unparsed(self):
        for rank, primer in enumerate(self.p3, start=1):
            print "Rank {0}\n-----------".format(rank)
            for attribute in primer.__dict__:
                print "{0}\t{1}".format(
                    attribute, primer.__dict__[attribute])
            print

    def fetch_top_primer(self):
        print "Top Forward: {0}\nTop Reverse: {1}\nTemplate: {2}\n".format(
            self.p3[0].forward_seq,
            self.p3[0].reverse_seq,
            self.full_seq)

    def run(self):
        '''Run the Illumina Primer creator'''

        print "\nYour run parameters:\
              \n--------------------\
              \nFull Region:\n{0}\
              \n\nRegion of Interest:\n{1}\
              \n\nPosible Insertion Lenght:{2}\
              \nRead Length:{3}\
              \nPaird-End Read Length:{4}\
              \n-----------------------\n".format(self.full_seq,
                                                  self.seq_to_include,
                                                  self.insertion,
                                                  self.region_length,
                                                  self.region_length / 2)
        for rank, primer in enumerate(self.p3, start=1):
            print "Rank: {0}\
                  \n---------".format(rank)
            d = self.get_primer_info(primer)
            for parameter in d:
                print "{0}: {1}".format(parameter, d[parameter])
            print

    def get_primer_info(self, primer):
        '''Where the magic happens, does the overlap conversions, 
        returns all the information in a dicitonary'''

        # First get all the forward stuff
        forward_primer = primer.forward_seq
        forward_gc = primer.forward_gc
        forward_tm = primer.forward_tm

        # Then get all the reverse stuff
        reverse_primer = primer.reverse_seq
        reverse_gc = primer.reverse_gc
        reverse_tm = primer.reverse_tm

        # The actual illumina read length
        illumina_read_length = self.region_length / 2

        # The possible insertions that can happen given area of interest
        possible_insertion = self.insertion

        # Magic
        # Get rev complement of returned primer
        rev_primer_rev_complement = str(Seq(
            reverse_primer, generic_dna).reverse_complement())

        # Index of region of interest
        _index_of_roi = self.full_seq.index(
            self.seq_to_include)

        # Index of forward primer
        _index_of_fwprimer = self.full_seq.index(
            forward_primer)

        '''The coverage is the index of the forward primer added 
        to a paird end illumina read length - 
        the index of the region of interest, draw it out!'''
        forward_coverage = (
            _index_of_fwprimer + illumina_read_length) - _index_of_roi

        # Index of the end of the region of interest
        _index_of_end_of_roi = _index_of_roi + len(self.seq_to_include)

        # Index of the end of the reverse primer
        _index_of_revprimer = self.full_seq.index(
            rev_primer_rev_complement) + len(reverse_primer)

        '''The reverse coverage is the illumina read length subtracted
        from the index of the reverse primer subtracted from the end 
        of the region of interest '''
        reverse_coverage = _index_of_end_of_roi - \
            (_index_of_revprimer - illumina_read_length)

        '''The overlap is the sequence of interest length plus the 
        possible insertions, subtracted from both of the coverages'''
        min_overlap = (
            forward_coverage + reverse_coverage) - len(self.seq_to_include)
        max_overlap = min_overlap + possible_insertion

        if min_overlap < 0:
            return False

        return_dict = OrderedDict()

        return_dict['Forward Primer'] = forward_primer
        return_dict['Forward Primer Length'] = len(forward_primer)
        return_dict['Forward GC Content'] = forward_gc
        return_dict['Forward Primer Tm'] = forward_tm
        return_dict['Reverse Primer'] = reverse_primer
        return_dict['Reverse Primer Length'] = len(reverse_primer)
        return_dict['Reverse GC Content'] = reverse_gc
        return_dict['Reverse Primer Tm'] = reverse_tm
        return_dict['Forward Region of Interest Coverage'] = forward_coverage
        return_dict['Reverse Region of Interest Coverage'] = reverse_coverage
        return_dict['Minimum Paired-End Overlap'] = min_overlap
        return_dict['Maximum Paired-End Overlap'] = max_overlap
        return return_dict

if __name__ == '__main__':
    IlluminaPrimer().run()
