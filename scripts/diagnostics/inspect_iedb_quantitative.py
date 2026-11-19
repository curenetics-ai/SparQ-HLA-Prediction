import requests
import json

url = "https://query-api.iedb.org/mhc_search"

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
        "reference_id"
    ),
    "limit": 5
}

try:
    response = requests.get(
        url,
        params=params,
        timeout=(15, 120)
    )

    print("Status:", response.status_code)
    response.raise_for_status()

    records = response.json()

    print("Records returned:", len(records))

    for record in records:
        print(json.dumps(record, indent=2))

except requests.exceptions.Timeout:
    print("IEDB request timed out. No dataset downloaded.")

except requests.exceptions.RequestException as error:
    print("Request failed:", error)