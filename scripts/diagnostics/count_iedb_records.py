
import requests

URL = "https://query-api.iedb.org/mhc_search"

headers = {
    "Range-Unit": "items",
    "Range": "0-0",
    "Prefer": "count=planned"
}

base_params = {
    "structure_type": "eq.Linear peptide",
    "mhc_class": "eq.I",
    "assay_names": "ilike.*binding*",
    "select": "elution_id"
}

for locus in ["A", "B", "C"]:

    params = base_params.copy()
    params["mhc_allele_name"] = f"like.HLA-{locus}*"

    try:
        response = requests.get(
            URL,
            params=params,
            headers=headers,
            timeout=(15, 90)
        )

        response.raise_for_status()

        count = response.headers.get("Content-Range", "Unavailable")

        print(f"HLA-{locus}: {count}")

    except requests.exceptions.RequestException as error:
        print(f"HLA-{locus}: Request failed - {error}")
