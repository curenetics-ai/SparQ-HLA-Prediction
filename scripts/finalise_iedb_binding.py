
import csv
from pathlib import Path
from collections import defaultdict, Counter

SOURCE = Path(
    "data/processed/binding/iedb_2chain_pair_status.csv"
)

OUTPUT = Path("data/processed/binding")
TABLES = Path("results/tables")

OUTPUT.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

stats = defaultdict(Counter)
training_rows = []
conflict_rows = []

with SOURCE.open(encoding="utf-8-sig") as file:

    reader = csv.DictReader(file)
    original_fields = reader.fieldnames

    for row in reader:

        allele = row["hla_allele"]
        status = row["classification"]

        stats[allele][status] += 1
        stats[allele]["total"] += 1

        if status == "POSITIVE":
            row["binary_label"] = 1
            training_rows.append(row)

        elif status == "NEGATIVE":
            row["binary_label"] = 0
            training_rows.append(row)

        elif status == "CONFLICTING":
            conflict_rows.append(row)

        else:
            raise ValueError(
                f"Unexpected classification: {status}"
            )


# Save unambiguous training candidates

training_path = OUTPUT / "iedb_2chain_training_candidates.csv"

with training_path.open(
    "w", newline="", encoding="utf-8-sig"
) as file:

    writer = csv.DictWriter(
        file,
        fieldnames=original_fields + ["binary_label"]
    )

    writer.writeheader()
    writer.writerows(training_rows)


# Save final allele-level binder summary

summary = []

for allele, counts in sorted(stats.items()):

    summary.append({
        "hla_allele": allele,
        "total_unique_pairs": counts["total"],
        "unambiguous_positive_binders": counts["POSITIVE"],
        "unambiguous_negative_pairs": counts["NEGATIVE"],
        "conflicting_pairs": counts["CONFLICTING"],
        "eligible_training_pairs": (
            counts["POSITIVE"] + counts["NEGATIVE"]
        )
    })

summary_path = TABLES / "final_binder_counts_by_allele.csv"

with summary_path.open(
    "w", newline="", encoding="utf-8-sig"
) as file:

    writer = csv.DictWriter(
        file,
        fieldnames=summary[0].keys()
    )

    writer.writeheader()
    writer.writerows(summary)


# Report

print("\n========== FINAL DATASET ==========")
print("Alleles:", len(summary))
print("Training candidates:", len(training_rows))
print("Positive binders:", sum(
    c["POSITIVE"] for c in stats.values()
))
print("Negative pairs:", sum(
    c["NEGATIVE"] for c in stats.values()
))
print("Conflicting pairs:", len(conflict_rows))

print("\nSaved files:")
print(training_path)
print(summary_path)
