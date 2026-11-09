"""Vendored from DMS_seq/src/script/reference_data.py (lab code), unchanged."""


#from src.script.read_bam import from_roman
import pandas as pd
from pandas import Series

def to_roman(number, skip_error=False):
    """
    Convert number to roman numeral
    """
    try:
        return {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI', 7: 'VII',
         8: 'VIII', 9: 'IX', 10: 'X', 11: 'XI', 12: 'XII', 13: 'XIII', 
         14: 'XIV', 15: 'XV', 16: 'XVI'}[number]
    except KeyError:
        if skip_error: return -1
        else: raise ValueError(f"Unsupported value: {number}")
    

def from_roman(roman, skip_error=False):
    """
    Convert Roman numeral to number
    """
    try:
        return {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, 
            "VI": 6, "VII": 7, "VIII": 8, "IX": 9, "X": 10, 
            "XI": 11, "XII": 12, "XIII": 13, "XIV": 14, 
            "XV": 15, "XVI": 16}[roman]
    except KeyError:
        if skip_error: return -1
        else: raise ValueError(f"Unsupported value: {roman}")

def parse_roman_chrom_series(data_chr):
    """Parse a Pandas Series which has chromosomes labeled with roman numerals
    into a list of integers. e.g. [chrI, chrIV] to [1, 4]"""
    chrom = data_chr.str.replace('chr', '').apply(lambda x: from_roman(x, True))
    return chrom


def read_sgd(filename):
    """Read SGD data file as tsv file from gff file with fasta data removed.
    This function returns all annotated elements in the genome."""

    data = pd.read_csv(filename, sep='\t', skiprows=19, 
                              names=["chr", "source", "cat", "start", "stop", ".", 
                              "strand", "", "desc"])
    data = data[data.columns[[0, 2, 3, 4, 6, 8]]]
    data.columns = ["chr", "cat", "start", "stop", "strand", "desc"]
    chrom = parse_roman_chrom_series(data.chr)
    data.chr = chrom

    return data


def extract_desc_val(data, key):
    """Extract the values in the description field of SGD"""
    def _extract_desc_val_row(orf_row, key):
        """Parse the description of a orf_row in sgd gff, extract `keys` and return 
        as a series. Can be used to return as a df when called with apply"""

        des_map = {}
        for entry in orf_row.desc.split(';'):
            k, val = tuple(entry.split('='))
            if k == key:
                return val

        return None

    # parse description column to extract relevant columns
    vals = data.apply(lambda row: _extract_desc_val_row(row, key=key),
        axis=1)

    return vals


def read_sgd_orfs(filename):
    """
    Read all ORFs from the SGD gff file
    """

    data = read_sgd(filename)
    orfs = data[(data['cat'] == 'gene') | (data['cat'] == 'rDNA')].copy()
    orfs['orf_name'] = None
    orfs['orf_class'] = None

    # chromosomal orfs
    orfs = orfs[orfs.chr > 0]

    orf_names = extract_desc_val(orfs, 'ID') 
    orf_classes = extract_desc_val(orfs, 'orf_classification')
    names = extract_desc_val(orfs, 'gene') 
    ontology = extract_desc_val(orfs, 'Ontology_term') 

    orfs['orf_name'] = orf_names
    orfs['orf_class'] = orf_classes

    # set name if it exists, or orf_name by default
    orfs['name'] = names
    orfs.loc[orfs.name.isna(), 'name'] = orfs[orfs.name.isna()]['orf_name']
    orfs['ontology'] = ontology

    orfs['length'] = orfs['stop'] - orfs['start'] + 1
    orfs = orfs.set_index('orf_name')[[
        'name','chr','start','stop','length','strand','orf_class', 'ontology']]

    return orfs


def read_sgd_orf_introns(filename):
    """Read all ORF intron data, useful for plotting introns within gene ORFS"""

    genes = read_sgd_orfs(filename)
    data = read_sgd(filename)

    keep_classes = ['intron', 'five_prime_UTR_intron', 'CDS']
    introns_CDSs = data[data.cat.isin(keep_classes)].copy()

    # rename 5' class to intron
    introns_CDSs.loc[(introns_CDSs['cat'] == 'five_prime_UTR_intron'), 'cat'] = 'intron'
    introns_CDSs['parent'] = extract_desc_val(introns_CDSs, 'Parent') 

    # relevant genes
    introns_CDSs = introns_CDSs[introns_CDSs.parent.isin(genes.index.values)]

    return introns_CDSs[['cat', 'start', 'stop', 'parent']].copy().reset_index(drop=True)


def read_park_TSS_PAS(TSS_filename, PAS_filename):
    """Read the Park 2014 TSS and PAS/TES dataset"""

    def _read_park_gff(filename, key):
        park_data = pd.read_csv(filename, sep='\t', skiprows=3)
        columns = park_data.columns[[0, 2, 3]] # relevant columns
        park_data = park_data[columns].copy()

        # cleanup
        park_data.columns = ['chr', 'orf_name', key]
        park_data.chr = parse_roman_chrom_series(park_data.chr)
        park_data.orf_name = park_data.orf_name.str.replace('_%s' % key, '')
        park_data[key] = park_data[key].astype(int)

        return park_data.set_index('orf_name')

    TSS = _read_park_gff(TSS_filename, 'TSS')
    PAS = _read_park_gff(PAS_filename, 'PAS')

    data = TSS[['TSS']].join(PAS[['PAS']])

    return data
