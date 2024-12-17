import os
from functools import lru_cache

from .helpers import csv_string_to_codons_dict, table_with_U_replaced_by_T
from .python_codon_tables_cocoputs import get_codons_table_cocoputs
from .python_codon_tables_kazusa import download_codons_table_kazusa


KNOWN_CUTS_SOURCES = {"kazusa", "cocoputs"}

_this_dir = os.path.dirname(os.path.realpath(__file__))
_tables_dir = os.path.join(_this_dir, "..", "codon_usage_data", "tables")

available_codon_tables_names = {
    source: [filename[:-4] for filename in os.listdir(os.path.join(_tables_dir, source))]
    for source in KNOWN_CUTS_SOURCES
}

available_codon_tables_shortnames = {
    source: {
        "_".join(table_name.split("_")[:-1]): table_name
        for table_name in available_codon_tables_names[source]
    }
    for source in KNOWN_CUTS_SOURCES
}


@lru_cache(maxsize=128)
def get_codons_table(table_name, replace_U_by_T=True, web_timeout=5, source="kazusa"):
    """Get data from one of this package's builtin codon usage tables.

    The ``table_name`` argument very flexible on purpose, it can be either an
    integer representing a taxonomic ID (which will be downloaded from
    the kazusa database), or a string "12245" representing a TaxID, or a string
    "e_coli_316407" referring to a builtin table of python_codon_optimization,
    or a short form "e_coli" which will be automatically extended to
    "e_coli_316407" (at your own risks).

    If a taxonomic ID is provided and no table with this taxID is present in
    the ``codon_usage_data/tables/`` folder, the table will be downloaded from
    either the http://www.kazusa.or.jp/codon or CocoPUTS
    (https://dnahive.fda.gov/dna.cgi?cmd=cuts_main) website. As these websites
    sometimes go down, the parameter ``web_timeout`` controls how long to wait
    before a Python exception is raised, informing the user that the service may
    be down.

    The ``replace_U_by_T`` argument will replace all codons names from UAA to
    TAA etc.

    The ``source`` argument can be used to specify the source of the codon tables.
    Either "kazusa" or "cocoputs" are supported (defaults to "kazusa").

    Returns a dict {"*": {'TAA': 0.64...}, 'K': {'AAA': 0.76...}, ...}

    
    """
    source = source.lower()
    if source not in KNOWN_CUTS_SOURCES:
        raise ValueError(f"Unknown source: {source}. Please select from {', '.join(KNOWN_CUTS_SOURCES)}")

    if replace_U_by_T:
        table = get_codons_table(table_name, replace_U_by_T=False, web_timeout=5, source=source)
        return table_with_U_replaced_by_T(table)
    if isinstance(table_name, int) or str.isdigit(table_name):
        if source == "kazusa":
            return download_codons_table_kazusa(taxid=table_name, timeout=web_timeout)
        elif source == "cocoputs":
            return get_codons_table_cocoputs(taxid=int(table_name))
    elif table_name in available_codon_tables_shortnames[source]:
        table_name = available_codon_tables_shortnames[source][table_name]
        with open(os.path.join(_tables_dir, source, table_name + ".csv"), "r") as f:
            return csv_string_to_codons_dict(f.read())
    else:
        raise ValueError(f"Unknown table_name: {table_name}")


def get_all_available_codons_tables(replace_U_by_T=True, source="kazusa"):
    """Get all data from all of this package's builtin codon usage tables."""
    return {
        table_name: get_codons_table(table_name, replace_U_by_T=replace_U_by_T, source=source)
        for table_name in available_codon_tables_names
    }
