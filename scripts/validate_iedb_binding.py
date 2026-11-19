
import json
import csv
from pathlib import Path
from collections import Counter

BASE = Path("data/raw/binding")
OUTPUT = Path("results/tables")

EXPECTED = {
    "HLA_A_2_chain": 22503,
    "HLA_A_class": 2652,
    "HLA_B_2_chain": 10565,
    "HLA_B_class": 851,
    "HLA_C_2_chain": 1804,
    "HLA_C_class": 20,
}

OUTPUT.mkdir(parents=True, exist_ok=True)

all_ids = set()
summary = []
total_records = 0
total_exports = 0
total_empty_exports = 0
errors = []

for partition, expected_count in EXPECTED.items():

    folder = BASE / partition
    files = sorted(folder.glob("page_*.jsonl"))

    count = 0
    peptides = set()
    peptide_alleles = set()
    alleles = set()

    outcomes = Counter()
    assay_types = Counter()
    units = Counter()

    empty_exports = 0
    export_count = 0

    locus = partition.split("_")[1]
    resolution = (
        "2 chain" if partition.endswith("2_chain")
        else "class"
    )

    for file_path in files:

        with file_path.open(
            "r", encoding="utf-8"
        ) as file:

            for line_number, line in enumerate(file, 1):

                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    errors.append(
                        f"Malformed JSON: {file_path}, "
                        f"line {line_number}"
                    )
                    continue

                count += 1
                total_records += 1

                record_id = record.get("elution_id")
                allele = record.get("mhc_allele_name")
                peptide = record.get("linear_sequence")

                if not record_id:
                    errors.append(
                        f"Missing elution_id: {file_path}"
                    )
                elif record_id in all_ids:
                    errors.append(
                        f"Duplicate elution_id: {record_id}"
                    )
                else:
                    all_ids.add(record_id)

                if record.get("mhc_class") != "I":
                    errors.append(
                        f"Incorrect MHC class: {record_id}"
                    )

                if not allele or not allele.startswith(
                    f"HLA-{locus}"
                ):
                    errors.append(
                        f"Incorrect locus: {record_id}"
                    )

                if record.get(
                    "mhc_allele_resolution"
                ) != resolution:
                    errors.append(
                        f"Incorrect resolution: {record_id}"
                    )

                if not peptide:
                    errors.append(
                        f"Missing peptide: {record_id}"
                    )
                else:
                    peptides.add(peptide)
                    peptide_alleles.add(
                        (peptide, allele)
                    )

                if allele:
                    alleles.add(allele)

                exports = record.get("mhc_export", [])

                if not exports:
                    empty_exports += 1

                for assay in exports:

                    export_count += 1
                    total_exports += 1

                    outcomes[
                        str(assay.get(
                            "assay__qualitative_measurement"
                        ))
                    ] += 1

                    assay_types[
                        str(assay.get(
                            "assay__response_measured"
                        ))
                    ] += 1

                    units[
                        str(assay.get("assay__units"))
                    ] += 1

    total_empty_exports += empty_exports

    if count != expected_count:
        errors.append(
            f"{partition}: expected {expected_count}, "
            f"found {count}"
        )

    print(f"\n--- {partition} ---")
    print("Records:", count)
    print("Unique alleles:", len(alleles))
    print("Unique peptides:", len(peptides))
    print("Unique peptide-allele pairs:", len(peptide_alleles))
    print("Nested assay records:", export_count)
    print("Empty export records:", empty_exports)
    print("Outcomes:", dict(outcomes))

    summary.append({
        "partition": partition,
        "expected_records": expected_count,
        "actual_records": count,
        "unique_alleles": len(alleles),
        "unique_peptides": len(peptides),
        "unique_peptide_allele_pairs": len(peptide_alleles),
        "nested_assay_records": export_count,
        "empty_export_records": empty_exports,
        "qualitative_outcomes": json.dumps(dict(outcomes)),
        "assay_types": json.dumps(dict(assay_types)),
        "measurement_units": json.dumps(dict(units)),
    })

summary_path = OUTPUT / "iedb_raw_validation.csv"

with summary_path.open(
    "w", newline="", encoding="utf-8"
) as file:

    writer = csv.DictWriter(
        file,
        fieldnames=summary[0].keys()
    )
    writer.writeheader()
    writer.writerows(summary)

print("\n========== FINAL VALIDATION ==========")
print("Total records:", total_records)
print("Unique record IDs:", len(all_ids))
print("Nested assay records:", total_exports)
print("Records with empty exports:", total_empty_exports)
print("Validation errors:", len(errors))

if errors:
    for error in errors[:30]:
        print("ERROR:", error)
else:
    print("All structural checks passed.")

print("Summary saved:", summary_path)
