#!/usr/bin/env python
"""Retrieve a Codon usage table from kazusa.or.jp, and store in a CSV file.

Usage:
------
To retrieve a codon table for one organism, given its TaxID, use:

> python codon_usage_retriever.py [TaxidNumber] [TargetFile.csv]

For instance:

> python codon_usage_retriever.py 316407 e_coli_codon_usage.csv

To retrieve codon tables from all organisms in ``organisms.csv`` at once, use:

> python codon_usage_retriever.py all



"""
import sys
import os
from python_codon_tables.python_codon_tables_kazusa import download_codons_table_kazusa
from python_codon_tables.python_codon_tables_cocoputs import fetch_codon_usage_from_cocoputs_api

_this_dir = os.path.dirname(os.path.realpath(__file__))

def _download_table_cocoputs(organism, taxid):
    target = os.path.join(_this_dir, "tables", "cocoputs", "%s_%s.csv" % (organism, taxid))
    result = fetch_codon_usage_from_cocoputs_api(taxid)
    result.to_csv(target)


def _download_table_kazusa(organism, taxid):
    target = os.path.join(_this_dir, "tables", "kazusa", "%s_%s.csv" % (organism, taxid))
    download_codons_table_kazusa(taxid=taxid, target_file=target)


def download_all_tables():
    with open(os.path.join(_this_dir, "organisms.csv"), "r") as f:
        for line in f.readlines()[1:]:
            organism, taxid = line.strip("\n").split(",")
            print("Retrieving %s (taxid %s)" % (organism, taxid))
            _download_table_kazusa(organism, taxid)
            _download_table_cocoputs(organism, taxid)


if __name__ == "__main__":
    print (" ".join(sys.argv))
    if sys.argv[1] == "all":
        download_all_tables()
    else:
        download_codons_table_kazusa(sys.argv[1], sys.argv[2])
        fetch_codon_usage_from_cocoputs_api(sys.argv[1]).to_csv(sys.argv[2])
