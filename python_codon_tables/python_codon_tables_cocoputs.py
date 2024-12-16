from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union
import requests

# Standard genetic code (Nuclear, standard)
# Codon to Amino Acid mapping (single letter)
CODON_TO_AA = {
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


@dataclass(frozen=True)
class CodonUsageResult:
    """Result of querying the Hive service for codon usage data."""

    taxid: int
    num_codons: int
    num_cds: int
    gc_percentage_overall: float
    gc_position_by_codon_position: tuple[float, float, float]
    codon_usage_table: dict[str, int]

    def group_codon_usage_table_by_amino_acid_and_normalize(self) -> dict[str, dict[str, float]]:
        """
        Return a dictionary mapping amino acid single-letter code to another
        dictionary of {codon: fraction}, where fraction is the fraction of that
        amino acid's codons that this particular codon constitutes.
        Stop codons will be included with amino acid '*'.
        """
        # Group codons by amino acid and sum their counts
        aa_grouped: dict[str, dict[str, int]] = {}
        for codon, count in self.codon_usage_table.items():
            if codon in CODON_TO_AA:
                aa = CODON_TO_AA[codon]
                if aa not in aa_grouped:
                    aa_grouped[aa] = {}
                aa_grouped[aa][codon] = count

        # Compute fractions
        aa_fractions: dict[str, dict[str, float]] = {}
        for aa, codons_dict in aa_grouped.items():
            total_count = sum(codons_dict.values())
            if total_count > 0:
                aa_fractions[aa] = {c: (count / total_count) for c, count in codons_dict.items()}
            else:
                # Edge case: if total_count == 0 (unlikely), just zero out.
                aa_fractions[aa] = {c: 0.0 for c in codons_dict}

        return aa_fractions

    def to_csv(self, csv_filepath: Union[str, Path]):
        """Convert the codon usage table to a CSV file

        Args:
            csv_filepath (Union[str, Path]): The path to the CSV file to output
        """
        if isinstance(csv_filepath, str):
            csv_filepath = Path(csv_filepath)

        with csv_filepath.open("w") as fh:
            fh.write("amino_acid,codon,relative_frequency\n")
            codon_usage_table_by_aa = self.group_codon_usage_table_by_amino_acid_and_normalize()

            for aa in sorted(codon_usage_table_by_aa.keys()):
                for codon in sorted(codon_usage_table_by_aa[aa].keys()):
                    fh.write(f"{aa},{codon},{codon_usage_table_by_aa[aa][codon]}\n")

    @classmethod
    def from_hive_response(cls, response_text: str) -> CodonUsageResult:
        """A factory method to create a CodonUsageResult object from the response text of a Hive query"""
        lines = response_text.strip().split("\n")

        metadata = {}
        codon_usage = {}

        for line in lines[1:]:  # skip the header line
            key, value = line.split(",", 1)
            key = key.strip().strip('"')
            value = value.strip().strip('"')

            if len(key) == 3 and all(n in "ACGT" for n in key.upper()):
                codon_usage[key.upper()] = int(value)
            else:
                metadata[key] = value

        taxid = int(metadata["taxid"])
        num_codons = int(metadata["#codon"])
        num_cds = int(metadata["#CDS"])
        gc_overall = float(metadata["GC%"])
        gc1 = float(metadata["GC1%"])
        gc2 = float(metadata["GC2%"])
        gc3 = float(metadata["GC3%"])

        return cls(
            taxid=taxid,
            num_codons=num_codons,
            num_cds=num_cds,
            gc_percentage_overall=gc_overall,
            gc_position_by_codon_position=(gc1, gc2, gc3),
            codon_usage_table=codon_usage,
        )


def fetch_codon_usage_from_hive(taxid: int) -> CodonUsageResult:
    """Fetch codon usage data from the FDA DNA Hive service.

    Args:
        taxid (int): The taxonomy ID of the organism to fetch

    Returns:
        CodonUsageResult: The result of the query
    """
    url = "https://dnahive.fda.gov/dna.cgi"
    params = {
        "cmd": "ionTaxidCollapseExt",
        "svcType": "svc-codon-usage",
        "objId": "537",
        "fileSource": "Refseq_species.tsv",
        "plen": "3",
        "taxid": str(taxid),
        "filterInColName": '["Organelle"]',
        "filterIn": '["genomic"]',
        "searchDeep": "true",
        "raw": "1",
    }

    response = requests.get(url, params=params)
    response.raise_for_status()
    return CodonUsageResult.from_hive_response(response.text)


def get_codons_table_from_hive(taxid: int) -> dict[str, dict[str, float]]:
    """Get the codon usage table from the FDA DNA Hive service and group by amino acid."""
    result = fetch_codon_usage_from_hive(taxid)
    return result.group_codon_usage_table_by_amino_acid_and_normalize()
