import requests
import json

URL = "https://query-api.iedb.org/mhc_search"

params = {
    "structure_type": "eq.Linear peptide",
    "mhc_class": "eq.I",
    "mhc_allele_name": "eq.HLA-A*02:01",
    "assay_names": "ilike.*binding*",
    "quantitative_measure": "not.is.null",

    "select": (
        "elution_id,"
        "linear_sequence,"
        "mhc_allele_name,"
        "mhc_class,"
        "assay_names,"
        "qualitative_measure,"
        "quantitative_measure,"
        "reference_id,"
        "mhc_export("
        "assay_id,"
        "assay__response_measured,"
        "assay__measurement_inequality,"
        "assay__quantitative_measurement,"
        "assay__units,"
        "assay__qualitative_measurement,"
        "reference__pmid"
        ")"
    ),

    "order": "elution_id.asc",
    "limit": 5
}

response = requests.get(URL, params=params, timeout=(15, 120))

print("Status:", response.status_code)
response.raise_for_status()

records = response.json()

print("Quantitative records:", len(records))

for record in records:
    print(json.dumps(record, indent=2))