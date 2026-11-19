
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
            "mhc_allele_name,"
            "mhc_class,"
            "qualitative_measure"
        ),
        "order": "elution_id.asc",
        "limit": 5
    }

    try:
        response = requests.get(
            URL, params=params, timeout=(15, 120)
        )
        response.raise_for_status()

        records = response.json()

        print(f"\n--- HLA-{locus} ---")
        print("Records:", len(records))

        for record in records:
            allele = record.get("mhc_allele_name")

            print(
                "Allele:", repr(allele),
                "| Class:", repr(record.get("mhc_class")),
                "| Peptide:", record.get("linear_sequence")
            )

            if not allele or not allele.startswith(f"HLA-{locus}"):
                print("CHECK: Unexpected allele naming")

            if record.get("mhc_class") != "I":
                print("CHECK: Unexpected MHC class")

    except requests.exceptions.RequestException as error:
        print(f"HLA-{locus} request failed:", error)

print("\nDiagnostic completed.")
