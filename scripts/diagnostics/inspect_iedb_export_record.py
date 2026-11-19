import requests
import json

url = "https://query-api.iedb.org/mhc_export"

params = {
    "assay_id": "eq.22666044",
    "select": (
        "assay_id,"
        "epitope__name,"
        "assay__response_measured,"
        "assay__method,"
        "assay__iri,"
        "assay__qualitative_measurement,"
        "assay__measurement_inequality,"
        "assay__quantitative_measurement,"
        "assay__units,"
        "mhc_restriction__name,"
        "mhc_restriction__class,"
        "reference__pmid"
    ),
    "limit": 1
}

response = requests.get(url, params=params, timeout=60)

print("Status:", response.status_code)
response.raise_for_status()

records = response.json()

print("Records returned:", len(records))

for record in records:
    print(json.dumps(record, indent=2))