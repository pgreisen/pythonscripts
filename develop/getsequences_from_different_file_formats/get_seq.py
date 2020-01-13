from Bio import SeqIO
def load_sequences_from_file(filename, nameToSeqOnly=True):
    """ Attempts to identify the type of sequence file, and return a dictionary
        with (id:sequence) pairs. Identifies file format by file extension.
        If nameToSeqOnly is False, maps to Bio.SeqRecord.SeqRecord objects instead """
    exts = {}; extToFlag = {}
    exts['fasta'] = ['fasta', 'fas', 'fa', 'faa', 'fna', 'frn']
    exts['fastq'] = ['fastq', 'fq']
    exts['genbank'] = ['genbank', 'gbk','gb']
    exts['embl'] = ['embl', 'em']
    exts['swiss'] = ['swiss', 'sw', 'swissprot']
    for flag in exts:
        for ext in exts[flag]:
            extToFlag[ext] = flag
    
    sequences = {}
    if os.path.isfile(filename):
        ext = filename.split('.')[-1].lower()
        if not ext in extToFlag: # not supported
            return {}
        records = SeqIO.parse(open(filename,'rU'), extToFlag[ext])
        for rec in records:
            if nameToSeqOnly:
                sequences[rec.id] = rec.seq
            else:
                sequences[rec.id] = rec
    return sequences
