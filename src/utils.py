import re

try:
    import pyhgvs as hgvs
    import pyhgvs.utils as hgvs_utils
    from pyfaidx import Fasta
except ImportError:
    Fasta = None


_COMP = dict(A='T', C='G', G='C', T='A', N='N',
             a='t', c='g', g='c', t='a', n='n')


def parse_hgvs(hgvs_name, fasta, genes):
    genome = Fasta(fasta, key_function=lambda x: 'chr{}'.format(x))

    with open(genes) as infile:
        transcripts = hgvs_utils.read_transcripts(infile)

    def get_transcript(name):
        return transcripts.get(name)

    return hgvs.parse_hgvs_name(hgvs_name, genome, get_transcript=get_transcript)


def getRevComp(seq):
    return ''.join(_COMP[base] for base in reversed(seq))


def getSequence(genome, chrom, start, end):
    return str(genome[str(chrom)][start - 1:end]).upper()


def parse_splice(cdsEffect, position, strand, fasta):
    genome = Fasta(fasta, key_function=lambda x: 'chr{}'.format(x))

    [chr, sPos] = position.split(':')
    startPos=int(sPos)
    if '>' in cdsEffect:
        mylist = cdsEffect.split('>')
        mylist[0] = re.sub(r'[A-Za-z]','',mylist[0])
        mytype = 'sub'
    elif '>' in cdsEffect:
        mylist = cdsEffect.split('>')
        mylist[0] = re.sub(r'[A-Za-z]','',mylist[0])
        mytype = 'sub'
    elif 'delins' in cdsEffect:
        mylist = cdsEffect.split('delins')
        mytype = 'delins'
    elif 'del' in cdsEffect:
        mylist = cdsEffect.split('del')
        mytype = 'del'
    elif 'ins' in cdsEffect:
        mylist = cdsEffect.split('ins')
        mytype = 'ins'
    elif 'dup' in cdsEffect:
        mylist = cdsEffect.split('dup')
        mytype = 'dup'
    else:
        raise ValueError('ERROR: not sure how to interpret [{}]')

    if strand == '-' and re.match(r'^[A-Za-z]+', mylist[1]):
        mylist[1] = getRevComp(mylist[1])
    myrange = mylist[0].split('_')
    if len(myrange) == 1:
        if mytype == 'sub':
            ref = getSequence(genome, chr, startPos, startPos)
            return (chr, startPos, ref, mylist[1])
        elif mytype == 'del':
            ref = getSequence(genome, chr, startPos, startPos+1)
            return (chr, startPos, ref, ref[0])
        elif mytype == 'dup':
            ref = getSequence(genome, chr, startPos, startPos+1)
            return (chr, startPos, ref[0], ref)
        else:
            raise ValueError('ERROR: not a range value and not a substitution or deletion [{}]'.format(cdsEffect))
    elif len(myrange) == 2:
        main = ['', '']
        sub = ['', '']
        for i in range(len(myrange)):
            if '+' in myrange[i]:
                [main[i], sub[i]] = re.split('\+',myrange[i])
            elif myrange[i].startswith('-'):
                main[i] = myrange[i]
                sub[i] = '0'
            elif '-' in myrange[i]:
                [main[i], sub[i]] = re.split('-',myrange[i])
                sub[i] = '-' + sub[i]
            else:
                main[i] = myrange[i]
                sub[i] = '0'

        if re.match(r'^[0-9]+', mylist[1]):
            seqlen = int(mylist[1])
        else:
            seqlen = len(mylist[1])

        if main[0].startswith('*') or main[1].startswith('*'):
            mylen = seqlen
        elif main[0] == main[1]:
            mylen = abs(int(sub[0]) - int(sub[1]))+1
        else:
            if sub[0] != '0' and sub[1] != '0':
                mylen = seqlen
            elif sub[0] == '0':
                mylen = abs(int(main[0]) - int(main[1]))+1
                mylen += abs(int(sub[1]))
            elif sub[1] == '0':
                mylen = abs(int(main[0]) - int(main[1]))+1
                mylen += abs(int(sub[0]))

        if mytype == 'ins':
            if mylen != 2:
                raise ValueError('ERROR: insertion but range is not 1 [{}]'.format(cdsEffect))
            elif not re.match(r'^[A-Za-z]+', mylist[1]):
                ref = getSequence(genome, chr, startPos, startPos)
                mylist[1] = ref + ('N' * int(mylist[1]))
            else:
                ref = getSequence(genome, chr, startPos, startPos)
                mylist[1] = ref + mylist[1]
        elif mytype == 'del':
        #    startPos = startPos - 1
        #    if mylen != seqlen:
        #        raise ValueError('ERROR: length of cds range does not match the given deleted sequence  [{}]'.format(cdsEffect))
            ref = getSequence(genome, chr, startPos, startPos+seqlen)
            mylist[1] = ref[0]
        elif mytype == 'dup':
            if mylen != seqlen:
                raise ValueError('ERROR: length of cds range does not match the given duplicated sequence  [{}]'.format(cdsEffect))
            ref = getSequence(genome, chr, startPos, startPos+mylen)
            mylist[1] = ref
            ref = ref[0]
        #elif main[0].startswith('*') or main[1].startswith('*'):
        #    raise ValueError('ERROR: insert+delete on 5UTR, unable to resolve sequence [{}]'.format(cdsEffect))
        else:
            ref = getSequence(genome, chr, startPos, startPos+mylen-1)

        return (chr, startPos, ref, mylist[1])
    else:
        raise ValueError('ERROR: syntax dont fit')
