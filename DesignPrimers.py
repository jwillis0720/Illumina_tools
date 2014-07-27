import argparse
import tempfile
import sys
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
        for primer in self.p3:
            self.get_primer_info(primer)

    def get_primer_info(self, primer):
        forward_read = primer.forward_seq
        forward_gc = primer.forward_gc
        forward_start = primer.forward_start
        forward_tm = primer.forward_tm
        forward_length = primer.forward_length

        reverse_read = primer.reverse_seq
        reverse_gc = primer.reverse_gc
        reverse_start = primer.reverse_start
        reverse_tm = primer.reverse_tm
        reverse_length = primer.reverse_length

        illumina_read_length = self.region_length / 2

        possible_insertion = self.insertion

        _index_of_roi = self.full_seq.index(self.seq_to_include)
        _index_of_fwprimer = self.full_seq.index(forward_read)
        _index_of_revprimer = self.full_seq.index(str(
            Seq(str(reverse_read), generic_dna).reverse_complement()))
        _index_of_end_of_roi = _index_of_roi + len(self.seq_to_include)

        forward_coverage = (
            _index_of_fwprimer + illumina_read_length) - _index_of_roi

        reverse_coverage = _index_of_end_of_roi - (
            _index_of_revprimer + len(reverse_read) - illumina_read_length)

        overlap = (
            forward_coverage + reverse_coverage) - (len(self.seq_to_include) + possible_insertion)

        print "Forward Primer:{0}\nForward Read Length:{1}\
               \nReverse Primer:{2}\nReverse Read Length:{3}\
               \nForward Coverage:{4}\
               \nReverse Coverage:{5}\
               \nMinimum Overlap: {6}".format(
            forward_read,
            illumina_read_length,
            reverse_read,
            illumina_read_length,
            forward_coverage,
            reverse_coverage,
            overlap)


if __name__ == '__main__':
    IlluminaPrimer().run()

        # aline = pw2.align.localxs(
        #    self.full_seq,
        #    p3[0].__dict__['forward_seq'],
        #    -10, -10)[0]
        # print aline[0] + "\n" + aline[1]
