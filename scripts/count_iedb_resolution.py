
import requests

URL = "https://query-api.iedb.org/mhc_search"

HEADERS = {
    "Range-Unit": "items",
    "Range": "0-0",
    "Prefer": "count=exact"
}

BASE = {
    "structure_type": "eq.Linear peptide",
    "mhc_class": "eq.I",
    "assay_names": "ilike.*binding*",
    "select": "elution_id"
}


def get_exact_count(params):
    response = requests.get(
        URL,
        params=params,
        headers=HEADERS,
        timeout=(15, 120)
    )
    response.raise_for_status()

    content_range = response.headers.get("Content-Range", "")
    total = content_range.split("/")[-1]

    if total == "*" or not total:
        raise ValueError(f"Exact count unavailable: {content_range}")

    return int(total)


for locus in ["A", "B", "C"]:

    params = BASE.copy()
    params["mhc_allele_name"] = f"like.HLA-{locus}*"

    total = get_exact_count(params)

    print(f"\n========== HLA-{locus} ==========")
    print("TOTAL:", total)

    categories = {}

    for resolution in ["2 chain", "1 chain", "class"]:

        query = params.copy()
        query["mhc_allele_resolution"] = f"eq.{resolution}"

        categories[resolution] = get_exact_count(query)

    query = params.copy()
    query["mhc_allele_resolution"] = "is.null"

    categories["NULL"] = get_exact_count(query)

    for category, count in categories.items():
        print(f"{category}: {count}")

    classified = sum(categories.values())

    print("Unaccounted records:", total - classified)

print("\nCount verification completed.")
