
import gzip
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

# ---------------- CONFIGURATION ----------------

URL = "https://query-api.iedb.org/mhc_search"

OUTPUT = Path("data/raw/processing/elution")
STATE_FILE = OUTPUT / "download_state.json"

PAGE_SIZE = 1000
CHUNK_SIZE = 10000
MAX_ATTEMPTS = 8

LOCI = ["A", "B", "C"]

RETRYABLE_STATUS = {408, 429, 500, 502, 503, 504}

FIELDS = (
    "elution_id,elution_iri,structure_id,structure_iri,"
    "linear_sequence,structure_type,linear_sequence_length,"
    "mhc_restriction,mhc_class,mhc_allele_name,"
    "mhc_allele_iri,mhc_allele_resolution,"
    "mhc_allele_evidence,qualitative_measure,"
    "assay_description,assay_names,assay_iris,"
    "reference_id,reference_iri,pubmed_id,"
    "reference_authors,reference_titles,reference_dates,"
    "journal_name,parent_source_antigen_iri,"
    "parent_source_antigen_name,"
    "parent_source_antigen_source_org_name,"
    "host_organism_name,source_organism_name"
)

OUTPUT.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({
    "User-Agent": "SparQ-HLA-Prediction/1.0"
})


# ---------------- NETWORK RETRIES ----------------

def request_json(params, headers, description):

    for attempt in range(1, MAX_ATTEMPTS + 1):

        try:
            response = session.get(
                URL,
                params=params,
                headers=headers,
                timeout=(30, 240)
            )

            if response.status_code in RETRYABLE_STATUS:

                delay = min(2 ** (attempt - 1), 60)

                print(
                    f"{description}: HTTP {response.status_code}. "
                    f"Retry {attempt}/{MAX_ATTEMPTS} "
                    f"in {delay}s."
                )

                response.close()
                time.sleep(delay)
                continue

            response.raise_for_status()

            records = response.json()

            if not isinstance(records, list):
                raise ValueError(
                    "Unexpected API response format."
                )

            return records, response.headers

        except (
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.ConnectionError,
            requests.exceptions.Timeout,
            requests.exceptions.ContentDecodingError,
            requests.exceptions.JSONDecodeError
        ) as error:

            if attempt == MAX_ATTEMPTS:
                raise

            delay = min(2 ** (attempt - 1), 60)

            print(
                f"{description}: Connection interrupted. "
                f"Retry {attempt}/{MAX_ATTEMPTS} "
                f"in {delay}s."
            )

            time.sleep(delay)

    raise RuntimeError(
        f"Maximum retries exceeded: {description}"
    )


# ---------------- CHECKPOINT MANAGEMENT ----------------

def save_state(state):

    temporary = STATE_FILE.with_suffix(".tmp")

    temporary.write_text(
        json.dumps(state, indent=2),
        encoding="utf-8"
    )

    temporary.replace(STATE_FILE)


if STATE_FILE.exists():

    state = json.loads(
        STATE_FILE.read_text(encoding="utf-8")
    )

else:

    state = {
        "next_offset": {},
        "completed": {}
    }

state.setdefault("next_offset", {})
state.setdefault("completed", {})


# ---------------- QUERY PARAMETERS ----------------

def base_params(locus):

    return {
        "structure_type": "eq.Linear peptide",
        "mhc_class": "eq.I",
        "mhc_allele_name": f"like.HLA-{locus}*",
        "mhc_allele_resolution": "eq.2 chain",
        "assay_names": "ilike.*mass spectrometry*",
        "select": FIELDS,
        "order": "elution_id.asc"
    }


# ---------------- EXACT RECORD COUNT ----------------

def exact_count(locus):

    headers = {
        "Range-Unit": "items",
        "Range": "0-0",
        "Prefer": "count=exact"
    }

    records, response_headers = request_json(
        base_params(locus),
        headers,
        f"HLA-{locus} count"
    )

    content_range = response_headers.get(
        "Content-Range", ""
    )

    if "/" not in content_range:
        raise ValueError(
            f"Missing exact count: {content_range}"
        )

    count = content_range.rsplit("/", 1)[1]

    if count == "*":
        raise ValueError(
            "IEDB did not return an exact count."
        )

    return int(count)


# ---------------- RECORD VALIDATION ----------------

def validate_record(record, locus):

    allele = record.get("mhc_allele_name", "")

    if record.get("mhc_class") != "I":
        raise ValueError("Unexpected MHC class.")

    if record.get("mhc_allele_resolution") != "2 chain":
        raise ValueError("Unexpected allele resolution.")

    if not allele.startswith(f"HLA-{locus}*"):
        raise ValueError(f"Unexpected allele: {allele}")

    if "mass spectrometry" not in str(
        record.get("assay_names", "")
    ).lower():
        raise ValueError("Unexpected assay type.")

    if not record.get("linear_sequence"):
        raise ValueError("Missing peptide sequence.")

    if not record.get("elution_id"):
        raise ValueError("Missing elution ID.")


# ---------------- DOWNLOAD EACH LOCUS ----------------

def download_locus(locus):

    folder = OUTPUT / f"HLA_{locus}_2_chain"
    folder.mkdir(parents=True, exist_ok=True)

    total = exact_count(locus)

    offset = int(
        state["next_offset"].get(locus, 0)
    )
   
    if offset > total or (
        offset != total and offset % CHUNK_SIZE != 0
    ):
        raise RuntimeError(
            f"Invalid checkpoint for HLA-{locus}: {offset}"
        )

    print(f"\nHLA-{locus}: {total:,} exact records")
    print(f"Resuming from: {offset:,}")

    while offset < total:

        chunk_end = min(
            offset + CHUNK_SIZE,
            total
        )

        chunk_number = offset // CHUNK_SIZE

        output_file = folder / (
            f"chunk_{chunk_number:05d}.jsonl.gz"
        )

        chunk_records = []
        cursor = offset

        # Retrieve one 10,000-record chunk.

        while cursor < chunk_end:

            end = min(
                cursor + PAGE_SIZE - 1,
                chunk_end - 1
            )

            headers = {
                "Range-Unit": "items",
                "Range": f"{cursor}-{end}",
                "Prefer": "count=exact"
            }

            params = base_params(locus)

            description = (
                f"HLA-{locus} records "
                f"{cursor:,}-{end:,}"
            )

            records, _ = request_json(
                params,
                headers,
                description
            )

            expected = end - cursor + 1

            if len(records) != expected:
                raise RuntimeError(
                    f"Expected {expected} records; "
                    f"received {len(records)}."
                )

            for record in records:
                validate_record(record, locus)

            chunk_records.extend(records)

            cursor += len(records)

            print(
                f"HLA-{locus}: "
                f"{cursor:,}/{total:,} retrieved"
            )

            time.sleep(0.1)

        if len(chunk_records) != chunk_end - offset:
            raise RuntimeError(
                "Incomplete chunk. Checkpoint not updated."
            )

        # Save completed chunk atomically.

        temporary_file = output_file.with_suffix(".tmp")

        with gzip.open(
            temporary_file,
            "wt",
            encoding="utf-8"
        ) as file:

            for record in chunk_records:

                file.write(
                    json.dumps(
                        record,
                        ensure_ascii=False
                    ) + "\n"
                )

        temporary_file.replace(output_file)

        # Update checkpoint only after successful save.

        offset = chunk_end

        state["next_offset"][locus] = offset
        state["completed"][locus] = offset == total

        save_state(state)

        print(
            f"Saved chunk {chunk_number:05d}: "
            f"{len(chunk_records):,} records"
        )

    print(f"Completed HLA-{locus}")


# ---------------- MAIN ----------------

def main():

    for locus in LOCI:
        download_locus(locus)

    manifest = {
        "project": "SparQ-HLA-Prediction",
        "source": "IEDB",
        "endpoint": URL,
        "extraction_date": datetime.now(
            timezone.utc
        ).isoformat(),
        "filters": {
            "structure_type": "Linear peptide",
            "mhc_class": "I",
            "loci": LOCI,
            "resolution": "2 chain",
            "assay": "Mass spectrometry"
        },
        "page_size": PAGE_SIZE,
        "chunk_size": CHUNK_SIZE,
        "download_state": state
    }

    manifest_path = OUTPUT / "manifest.json"

    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8"
    )

    print("\n========== EXTRACTION SUMMARY ==========")

    for locus in LOCI:
        print(
            f"HLA-{locus}: "
            f"{state['next_offset'].get(locus, 0):,}"
        )

    print("Manifest:", manifest_path)


if __name__ == "__main__":
    main()
