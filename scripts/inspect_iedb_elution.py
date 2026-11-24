
import requests
import json
from pathlib import Path

URL = "https://query-api.iedb.org/mhc_search"

OUT = Path("data/raw/processing/elution")
OUT.mkdir(parents=True, exist_ok=True)

for locus in ["A", "B", "C"]:

    params = {
        "structure_type": "eq.Linear peptide",
        "mhc_class": "eq.I",
        "mhc_allele_name": f"like.HLA-{locus}*",
        "mhc_allele_resolution": "eq.2 chain",
        "assay_names": "ilike.*mass spectrometry*",
        "select": (
            "elution_id,"
            "linear_sequence,"
            "mhc_class,"
            "mhc_allele_name,"
            "mhc_allele_resolution,"
            "assay_names,"
            "qualitative_measure,"
            "reference_id,"
            "reference_dates,"
            "parent_source_antigen_name,"
            "parent_source_antigen_iri"
        ),
        "order": "elution_id.asc",
        "limit": 5
    }

    headers = {
        "Range-Unit": "items",
        "Range": "0-4",
        "Prefer": "count=exact"
    }

    response = requests.get(
        URL,
        params=params,
        headers=headers,
        timeout=(15, 120)
    )

    response.raise_for_status()

    records = response.json()

    count = response.headers.get(
        "Content-Range", "Unavailable"
    )

    print(f"\nHLA-{locus}")
    print("Exact count:", count)
    print("Sample records:", len(records))

    output_file = OUT / f"pilot_HLA_{locus}.json"

    output_file.write_text(
        json.dumps(records, indent=2),
        encoding="utf-8"
    )

    for record in records[:2]:
        print(
            record["linear_sequence"],
            record["mhc_allele_name"],
            record["assay_names"]
        )

print("\nElution query inspection complete.")
