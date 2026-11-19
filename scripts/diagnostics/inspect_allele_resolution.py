
import requests

URL = "https://query-api.iedb.org/mhc_search"

for locus in ["A", "B", "C"]:

    params = {
        "structure_type": "eq.Linear peptide",
        "mhc_class": "eq.I",
        "mhc_allele_name": f"like.HLA-{locus}*",
        "assay_names": "ilike.*binding*",
        "select": (
            "elution_id,"
            "linear_sequence,"
            "mhc_restriction,"
            "mhc_allele_name,"
            "mhc_allele_iri,"
            "mhc_allele_resolution,"
            "mhc_allele_evidence,"
            "mhc_class"
        ),
        "order": "elution_id.asc",
        "limit": 10
    }

    response = requests.get(URL, params=params, timeout=(15, 120))
    response.raise_for_status()

    records = response.json()

    print(f"\n========== HLA-{locus} ==========")

    for record in records:
        print({
            "allele": record["mhc_allele_name"],
            "restriction": record["mhc_restriction"],
            "resolution": record["mhc_allele_resolution"],
            "evidence": record["mhc_allele_evidence"],
            "allele_iri": record["mhc_allele_iri"],
        })

print("\nResolution inspection completed.")
