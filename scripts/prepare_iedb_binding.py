
import csv
import json
from pathlib import Path
from collections import defaultdict, Counter

RAW = Path("data/raw/binding")
OUT = Path("data/processed/binding")
TABLES = Path("results/tables")

OUT.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

PARTITIONS = {
    "HLA_A_2_chain": 22503,
    "HLA_A_class": 2652,
    "HLA_B_2_chain": 10565,
    "HLA_B_class": 851,
    "HLA_C_2_chain": 1804,
    "HLA_C_class": 20,
}

POSITIVE = {
    "Positive",
    "Positive-Low",
    "Positive-Intermediate",
    "Positive-High",
}

observations = []
allele_stats = defaultdict(Counter)
pair_outcomes = defaultdict(set)

seen_ids = set()

for partition, expected in PARTITIONS.items():

    folder = RAW / partition
    files = sorted(folder.glob("page_*.jsonl"))

    count = 0

    locus = partition.split("_")[1]
    resolution = (
        "2 chain" if partition.endswith("2_chain")
        else "class"
    )

    for file_path in files:

        with file_path.open(
            encoding="utf-8"
        ) as file:

            for line in file:

                record = json.loads(line)

                count += 1

                record_id = record["elution_id"]

                if record_id in seen_ids:
                    raise ValueError(
                        f"Duplicate record ID: {record_id}"
                    )

                seen_ids.add(record_id)

                allele = record["mhc_allele_name"]
                peptide = record["linear_sequence"]

                exports = record["mhc_export"]

                if not exports:
                    raise ValueError(
                        f"Missing assay export: {record_id}"
                    )

                for assay in exports:

                    outcome = assay.get(
                        "assay__qualitative_measurement"
                    )

                    if outcome in POSITIVE:
                        status = "Positive"
                    elif outcome == "Negative":
                        status = "Negative"
                    else:
                        status = "Other/Unclassified"

                    row = {
                        "elution_id": record_id,
                        "assay_id": assay.get("assay_id"),
                        "peptide": peptide,
                        "peptide_length": len(peptide),
                        "hla_allele": allele,
                        "hla_locus": locus,
                        "mhc_class": record["mhc_class"],
                        "allele_resolution": resolution,
                        "qualitative_outcome": outcome,
                        "binary_outcome": status,
                        "assay_response": assay.get(
                            "assay__response_measured"
                        ),
                        "assay_method": assay.get(
                            "assay__method"
                        ),
                        "quantitative_measurement": assay.get(
                            "assay__quantitative_measurement"
                        ),
                        "measurement_units": assay.get(
                            "assay__units"
                        ),
                        "measurement_inequality": assay.get(
                            "assay__measurement_inequality"
                        ),
                        "reference_id": record.get(
                            "reference_id"
                        ),
                        "pmid": assay.get("reference__pmid"),
                    }

                    observations.append(row)

                    key = (resolution, allele)

                    allele_stats[key]["total"] += 1
                    allele_stats[key][status] += 1

                    # Unique peptide-allele observations
                    pair_key = (resolution, allele, peptide)
                    pair_outcomes[pair_key].add(status)

    if count != expected:
        raise ValueError(
            f"{partition}: expected {expected}, got {count}"
        )

    print(f"{partition}: {count:,} validated")


# Save the complete structured observations

observation_path = OUT / "iedb_binding_observations.csv"

fields = list(observations[0].keys())

with observation_path.open(
    "w", newline="", encoding="utf-8-sig"
) as file:

    writer = csv.DictWriter(file, fieldnames=fields)
    writer.writeheader()
    writer.writerows(observations)


# Aggregate unique peptide-allele outcomes

pair_stats = defaultdict(Counter)

for (resolution, allele, peptide), outcomes in pair_outcomes.items():

    stats = pair_stats[(resolution, allele)]

    stats["unique_pairs"] += 1

    if "Positive" in outcomes:
        stats["positive_pairs"] += 1

    if "Negative" in outcomes:
        stats["negative_pairs"] += 1

    if "Positive" in outcomes and "Negative" in outcomes:
        stats["conflicting_pairs"] += 1


# Create allele-level summary

summary = []

for (resolution, allele), counts in sorted(allele_stats.items()):

    pair = pair_stats[(resolution, allele)]

    summary.append({
        "allele": allele,
        "resolution": resolution,
        "total_assay_records": counts["total"],
        "positive_assay_records": counts["Positive"],
        "negative_assay_records": counts["Negative"],
        "other_records": counts["Other/Unclassified"],
        "unique_peptide_allele_pairs": pair["unique_pairs"],
        "unique_positive_peptide_pairs": pair["positive_pairs"],
        "unique_negative_peptide_pairs": pair["negative_pairs"],
        "conflicting_peptide_pairs": pair["conflicting_pairs"],
    })

summary_path = TABLES / "measured_binder_counts_by_allele.csv"

with summary_path.open(
    "w", newline="", encoding="utf-8-sig"
) as file:

    writer = csv.DictWriter(
        file,
        fieldnames=summary[0].keys()
    )
    writer.writeheader()
    writer.writerows(summary)


print("\n========== PROCESSING SUMMARY ==========")
print("Observation rows:", len(observations))
print("Unique elution IDs:", len(seen_ids))
print("Allele summary rows:", len(summary))
print("Primary two-chain allele groups:",
      sum(1 for r in summary if r["resolution"] == "2 chain"))
print("Class-level groups:",
      sum(1 for r in summary if r["resolution"] == "class"))

print("\nSaved:", observation_path)
print("Saved:", summary_path)
