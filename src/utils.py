try:
    import pyhgvs as hgvs
    import pyhgvs.utils as hgvs_utils
    from pyfaidx import Fasta
except ImportError:
    Fasta = None


def parse_hgvs(hgvs_name, fasta, genes):
    genome = Fasta(fasta, key_function = lambda x: 'chr{}'.format(x))

    with open(genes) as infile:
        transcripts = hgvs_utils.read_transcripts(infile)

    def get_transcript(name):
        return transcripts.get(name)

    return hgvs.parse_hgvs_name(hgvs_name, genome, get_transcript=get_transcript)
