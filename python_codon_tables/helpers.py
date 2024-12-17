# Standard genetic code (Nuclear, standard)
# Codon to Amino Acid mapping (single letter)
CODON_TO_AA_DNA = {
    # U -> T in DNA
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "TGC": "C",
    "TGT": "C",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "TTC": "F",
    "TTT": "F",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "CAC": "H",
    "CAT": "H",
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "AAA": "K",
    "AAG": "K",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "TTA": "L",
    "TTG": "L",
    "ATG": "M",
    "AAC": "N",
    "AAT": "N",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAA": "Q",
    "CAG": "Q",
    "AGA": "R",
    "AGG": "R",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "AGC": "S",
    "AGT": "S",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "TGG": "W",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",  # Stop codons
}


CODON_TO_AA_RNA = {codon.replace("T", "U"): amino_acid for codon, amino_acid in CODON_TO_AA_DNA.items()}


def csv_string_to_codons_dict(csv_string):
    """Transform a CSV string of a codon table to a dict."""
    result = {}
    for line in csv_string.split("\n")[1:]:
        aa, codon, freq = line.split(",")
        if aa not in result:
            result[aa] = {}
        result[aa][codon] = float(freq)
    return result


def table_with_U_replaced_by_T(table):
    return {
        aa: {codon.replace("U", "T"): freq for codon, freq in aa_data.items()}
        for aa, aa_data in table.items()
    }