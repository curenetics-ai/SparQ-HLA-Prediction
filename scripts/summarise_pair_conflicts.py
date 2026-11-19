
import json
import csv
from pathlib import Path
from collections import defaultdict, Counter

RAW = Path("data/raw/binding")
OUT = Path("data/processed/binding")
TABLES = Path("results/tables")

OUT.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

PARTITIONS = [
    "HLA_A_2_chain",
    "HLA_B_2_chain",
    "HLA_C_2_chain",
]

POSITIVE = {
    "Positive",
    "Positive-Low",
    "Positive-Intermediate",
    "Positive-High",
}

pairs = defaultdict(list)
seen_assay_ids = set()
duplicate_assay_ids = set()

for partition in PARTITIONS:

    for path in sorted((RAW / partition).glob("page_*.jsonl")):

        with path.open(encoding="utf-8") as file:

            for line in file:

                record = json.loads(line)

                allele = record["mhc_allele_name"]
                peptide = record["linear_sequence"]

                for assay in record.get("mhc_export", []):

                    assay_id = assay.get("assay_id")

                    if assay_id is not None:
                        if assay_id in seen_assay_ids:
                            duplicate_assay_ids.add(assay_id)
                        seen_assay_ids.add(assay_id)

                    outcome = assay.get(
                        "assay__qualitative_measurement"
                    )

                    if outcome in POSITIVE:
                        status = "Positive"
                    elif outcome == "Negative":
                        status = "Negative"
                    else:
                        status = "Other"

                    pairs[(allele, peptide)].append({
                        "assay_id": assay_id,
                        "outcome": outcome,
                        "status": status,
                        "method": assay.get("assay__method"),
                        "pmid": assay.get("reference__pmid"),
                    })


# Classify every unique peptide-allele pair

pair_rows = []
allele_summary = defaultdict(Counter)

for (allele, peptide), evidence in sorted(pairs.items()):

    positive = sum(
        item["status"] == "Positive" for item in evidence
    )

    negative = sum(
        item["status"] == "Negative" for item in evidence
    )

    other = sum(
        item["status"] == "Other" for item in evidence
    )

    if positive and negative:
        classification = "CONFLICTING"
    elif positive:
        classification = "POSITIVE"
    elif negative:
        classification = "NEGATIVE"
    else:
        classification = "UNCLASSIFIED"

    row = {
        "hla_allele": allele,
        "peptide": peptide,
        "peptide_length": len(peptide),
        "positive_assays": positive,
        "negative_assays": negative,
        "other_assays": other,
        "classification": classification,
        "assay_ids": json.dumps([
            item["assay_id"] for item in evidence
        ]),
        "pmids": json.dumps([
            item["pmid"] for item in evidence if item["pmid"]
        ]),
        "evidence": json.dumps(evidence),
    }

    pair_rows.append(row)

    allele_summary[allele][classification] += 1
    allele_summary[allele]["total_pairs"] += 1


# Save complete pair-level dataset

pair_fields = list(pair_rows[0].keys())

with (OUT / "iedb_2chain_pair_status.csv").open(
    "w", newline="", encoding="utf-8-sig"
) as file:

    writer = csv.DictWriter(file, fieldnames=pair_fields)
    writer.writeheader()
    writer.writerows(pair_rows)


# Save conflicting pairs separately

conflicts = [
    row for row in pair_rows
    if row["classification"] == "CONFLICTING"
]

with (TABLES / "iedb_2chain_conflicts.csv").open(
    "w", newline="", encoding="utf-8-sig"
) as file:

    writer = csv.DictWriter(file, fieldnames=pair_fields)
    writer.writeheader()
    writer.writerows(conflicts)


# Save revised allele summary

summary_rows = []

for allele, counts in sorted(allele_summary.items()):

    summary_rows.append({
        "hla_allele": allele,
        **dict(counts),
    })

summary_fields = list(summary_rows[0].keys())

with (TABLES / "iedb_2chain_pair_summary_by_allele.csv").open(
    "w", newline="", encoding="utf-8-sig"
) as file:

    writer = csv.DictWriter(file, fieldnames=summary_fields)
    writer.writeheader()
    writer.writerows(summary_rows)


# Final report

totals = Counter()

for row in pair_rows:
    totals[row["classification"]] += 1

print("\n========== PAIR CLASSIFICATION ==========")
print("Total unique peptide-allele pairs:", len(pair_rows))
print("Positive:", totals["POSITIVE"])
print("Negative:", totals["NEGATIVE"])
print("Conflicting:", totals["CONFLICTING"])
print("Unclassified:", totals["UNCLASSIFIED"])

print("\nUnique assay IDs:", len(seen_assay_ids))
print("Duplicate assay ID occurrences:", len(duplicate_assay_ids))

print("\nFiles saved:")
print(OUT / "iedb_2chain_pair_status.csv")
print(TABLES / "iedb_2chain_conflicts.csv")
print(TABLES / "iedb_2chain_pair_summary_by_allele.csv")
