import re
import sys

if sys.version_info[0] == 3:
    import urllib.request

    urlopen = urllib.request.urlopen
else:
    import urllib2

    urlopen = urllib2.urlopen

from functools import lru_cache

from .helpers import csv_string_to_codons_dict


@lru_cache(maxsize=128)
def download_codons_table_kazusa(taxid=316407, target_file=None, timeout=5):
    """Get all data from all of this package's builtin codon usage tables."""
    _kazusa_url = (
        "http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi" "?aa=1&style=N&species=%s"
    )
    _codon_regexpr = r"([ATGCU]{3}) ([A-Z]|\*) (\d.\d+)"
    url = _kazusa_url % taxid
    try:
        web_handle = urlopen(url, timeout=timeout)
    except Exception as err:
        if "timed out" in str(err):
            raise RuntimeError(
                (
                    "connexion to %s timed out after %d seconds. Maybe "
                    "their website is down?"
                )
                % (url, timeout)
            )
        else:
            raise err

    html_content = web_handle.read().decode().replace("\n", " ")
    if "<title>not found</title>" in html_content.lower():
        raise RuntimeError(
            "Codon usage table for taxonomy ID '%s' not found:" " %s" % (taxid, url)
        )
    csv_data = "\n".join(
        ["amino_acid,codon,relative_frequency"]
        + sorted(
            [
                "%s,%s,%s" % (aa, codon, usage)
                for codon, aa, usage in re.findall(_codon_regexpr, html_content)
            ]
        )
    )
    if target_file is not None:
        with open(target_file, "w+") as f:
            f.write(csv_data)
    else:
        return csv_string_to_codons_dict(csv_data)