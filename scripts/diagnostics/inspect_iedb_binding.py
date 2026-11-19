
import requests
import json

URL = "https://query-api.iedb.org/mhc_export"

BASE_PARAMS = {
    "assay__response_measured": "ilike.*binding*",
    "mhc_restriction__class": "eq.I",
    "mhc_restriction__name": "ilike.HLA-*",
    "select": (
        "assay_id,"
        "epitope__name,"
        "assay__response_measured,"
        "assay__method,"
        "assay__qualitative_measurement,"
        "assay__measurement_inequality,"
        "assay__quantitative_measurement,"
        "assay__units,"
        "mhc_restriction__name,"
        "mhc_restriction__class,"
        "reference__pmid"
    ),
    "order": "assay_id.desc",
    "limit": 5
}


def inspect_records(label, extra_filters=None):
    params = BASE_PARAMS.copy()

    if extra_filters:
        params.update(extra_filters)

    response = requests.get(URL, params=params, timeout=60)

    print(f"\n--- {label} ---")
    print("Status:", response.status_code)

    response.raise_for_status()

    records = response.json()
    print("Records:", len(records))

    for record in records:
        print(json.dumps(record, indent=2))

    return records


# Test 1: General Class I binding records
inspect_records("HLA Class I binding records")

# Test 2: Class I binding records with measured values
inspect_records(
    "Quantitative Class I binding records",
    {"assay__quantitative_measurement": "not.is.null"}
)
