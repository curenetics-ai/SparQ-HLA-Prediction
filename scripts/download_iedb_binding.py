
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# ---------------- CONFIGURATION ----------------

URL = "https://query-api.iedb.org/mhc_search"

OUTPUT_DIR = Path("data/raw/binding")

PAGE_SIZE = 500

LOCI = ["A", "B", "C"]

RESOLUTIONS = ["2 chain", "class"]


# ---------------- HTTP SESSION ----------------

retry_strategy = Retry(
    total=5,
    connect=5,
    read=5,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=frozenset(["GET"]),
    respect_retry_after_header=True
)

session = requests.Session()

session.mount(
    "https://",
    HTTPAdapter(max_retries=retry_strategy)
)

session.headers.update({
    "User-Agent": "SparQ-HLA-Prediction/1.0"
})


# ---------------- QUERY CONFIGURATION ----------------

def make_params(locus, resolution, select):
    return {
        "structure_type": "eq.Linear peptide",
        "mhc_class": "eq.I",
        "mhc_allele_name": f"like.HLA-{locus}*",
        "mhc_allele_resolution": f"eq.{resolution}",
        "assay_names": "ilike.*binding*",
        "select": select,
        "order": "elution_id.asc"
    }


def get_exact_count(params):

    headers = {
        "Range-Unit": "items",
        "Range": "0-0",
        "Prefer": "count=exact"
    }

    response = session.get(
        URL,
        params=params,
        headers=headers,
        timeout=(20, 120)
    )

    response.raise_for_status()

    content_range = response.headers.get(
        "Content-Range", ""
    )

    if "/" not in content_range:
        raise ValueError(
            f"Missing Content-Range: {content_range}"
        )

    total = content_range.rsplit("/", 1)[1]

    if total == "*":
        raise ValueError("Exact count was not returned.")

    return int(total)


# ---------------- PAGE VALIDATION ----------------

def valid_existing_page(path, expected_rows):

    if not path.exists():
        return False

    try:
        with path.open("r", encoding="utf-8") as file:
            records = [
                json.loads(line)
                for line in file
                if line.strip()
            ]

        return len(records) == expected_rows

    except (OSError, json.JSONDecodeError):
        return False


# ---------------- DATA EXTRACTION ----------------

def download_partition(locus, resolution):

    partition = f"HLA_{locus}_{resolution.replace(' ', '_')}"

    output = OUTPUT_DIR / partition

    output.mkdir(parents=True, exist_ok=True)

    count_params = make_params(
        locus,
        resolution,
        "elution_id"
    )

    total = get_exact_count(count_params)

    print(f"\n{partition}: {total:,} records")

    downloaded = 0
    export_missing = 0

    for offset in range(0, total, PAGE_SIZE):

        expected_rows = min(
            PAGE_SIZE,
            total - offset
        )

        page_number = offset // PAGE_SIZE

        page_path = output / (
            f"page_{page_number:05d}.jsonl"
        )

        if valid_existing_page(page_path, expected_rows):

            downloaded += expected_rows

            print(
                f"Existing page verified: "
                f"{downloaded:,}/{total:,}"
            )

            continue

        params = make_params(
            locus,
            resolution,
            "*,mhc_export(*)"
        )

        headers = {
            "Range-Unit": "items",
            "Range": (
                f"{offset}-{offset + PAGE_SIZE - 1}"
            ),
            "Prefer": "count=exact"
        }

        response = session.get(
            URL,
            params=params,
            headers=headers,
            timeout=(20, 120)
        )

        response.raise_for_status()

        records = response.json()

        if len(records) != expected_rows:
            raise RuntimeError(
                f"Unexpected page size at offset {offset}: "
                f"expected {expected_rows}, "
                f"received {len(records)}"
            )

        for record in records:

            allele = record.get("mhc_allele_name", "")

            if record.get("mhc_class") != "I":
                raise ValueError("Unexpected MHC class.")

            if not allele.startswith(f"HLA-{locus}"):
                raise ValueError(f"Unexpected allele: {allele}")

            if record.get("mhc_allele_resolution") != resolution:
                raise ValueError("Unexpected allele resolution.")

            if "mhc_export" not in record:
                raise ValueError("Missing export metadata.")

            if not record["mhc_export"]:
                export_missing += 1

        temp_path = page_path.with_suffix(".tmp")

        with temp_path.open("w", encoding="utf-8") as file:

            for record in records:

                file.write(
                    json.dumps(record, ensure_ascii=False)
                    + "\n"
                )

        temp_path.replace(page_path)

        downloaded += len(records)

        print(
            f"Downloaded: {downloaded:,}/{total:,}"
        )

        time.sleep(0.2)

    print(f"Completed: {partition}")

    return {
        "partition": partition,
        "exact_count": total,
        "downloaded": downloaded,
        "records_without_export": export_missing
    }


# ---------------- MAIN EXECUTION ----------------

def main():

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    results = []

    for locus in LOCI:

        for resolution in RESOLUTIONS:

            result = download_partition(
                locus,
                resolution
            )

            results.append(result)

    manifest = {
        "project": "SparQ-HLA-Prediction",
        "source": "IEDB IQ-API",
        "endpoint": URL,
        "extraction_date": datetime.now(
            timezone.utc
        ).isoformat(),
        "filters": {
            "structure_type": "Linear peptide",
            "mhc_class": "I",
            "loci": LOCI,
            "assay_names": "contains binding",
            "resolutions": RESOLUTIONS
        },
        "page_size": PAGE_SIZE,
        "partitions": results,
        "total_records": sum(
            item["downloaded"] for item in results
        )
    }

    manifest_path = OUTPUT_DIR / "manifest.json"

    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8"
    )

    print("\n========== EXTRACTION SUMMARY ==========")

    print(
        "Total downloaded:",
        f"{manifest['total_records']:,}"
    )

    print("Manifest:", manifest_path)


if __name__ == "__main__":
    main()
