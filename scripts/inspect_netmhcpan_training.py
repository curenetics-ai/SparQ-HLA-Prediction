
import csv
import re
from pathlib import Path
from collections import Counter

ROOT = Path("data/reference/netmhcpan41")

IEDB_FILE = Path(
    "data/processed/binding/iedb_2chain_pair_status.csv"
)

OUTPUT = Path("results/tables")
OUTPUT.mkdir(parents=True, exist_ok=True)


def normalise_allele(value):
    """Standardise HLA-A02:01 and HLA-A*02:01."""

    match = re.fullmatch(
        r"HLA-([ABC])\*?(\d{2}):?(\d{2})",
        value.strip(),
        flags=re.IGNORECASE
    )

    if not match:
        return None

    locus, first, second = match.groups()

    return f"HLA-{locus.upper()}*{first}:{second}"


# Load IEDB high-resolution peptide-allele pairs.

iedb_pairs = set()
iedb_sequences = set()

with IEDB_FILE.open(encoding="utf-8-sig") as file:
    for row in csv.DictReader(file):

        allele = normalise_allele(row["hla_allele"])
        peptide = row["peptide"].strip().upper()

        if allele and peptide:
            iedb_pairs.add((peptide, allele))
            iedb_sequences.add(peptide)


# Read multi-allelic cell-line mappings.

allele_files = list(ROOT.rglob("allelelist*"))

if not allele_files:
    raise FileNotFoundError("NetMHCpan allelelist not found.")

cell_mapping = {}

with allele_files[0].open(encoding="utf-8") as file:

    for line in file:

        parts = line.strip().split(maxsplit=1)

        if len(parts) != 2:
            continue

        cell_line, allele_string = parts

        mapped = {
            normalise_allele(item)
            for item in allele_string.split(",")
        }

        mapped.discard(None)
        cell_mapping[cell_line] = mapped


print("IEDB unique peptide-allele pairs:", len(iedb_pairs))
print("Cell-line mappings loaded:", len(cell_mapping))


inventory = []
overlap_summary = []
all_exact_overlaps = []
all_possible_overlaps = []

for dataset in ["ba", "el"]:

    direct_pairs = set()
    possible_pairs = set()
    reference_sequences = set()

    print(f"\n========== {dataset.upper()} DATA ==========")

    for fold in range(5):

        filename = f"c{fold:03d}_{dataset}"
        matches = list(ROOT.rglob(filename))

        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected one {filename}, found {len(matches)}"
            )

        file_path = matches[0]

        counts = Counter()
        file_pairs = set()

        with file_path.open(
            encoding="utf-8-sig"
        ) as file:

            for line in file:

                parts = line.split()

                if len(parts) < 3:
                    counts["malformed"] += 1
                    continue

                peptide = parts[0].upper()
                source = parts[2]

                counts["rows"] += 1
                reference_sequences.add(peptide)

                allele = normalise_allele(source)

                if allele:

                    counts["direct_hla_rows"] += 1
                    counts["hla_a_b_c_rows"] += 1

                    pair = (peptide, allele)

                    direct_pairs.add(pair)
                    file_pairs.add(pair)

                elif source in cell_mapping:

                    counts["cell_line_rows"] += 1

                    for allele in cell_mapping[source]:

                        counts["mapped_hla_associations"] += 1

                        pair = (peptide, allele)
                        possible_pairs.add(pair)

                elif source.upper().startswith("HLA-"):

                    counts["other_hla_designations"] += 1

                else:

                    counts["other_identifiers"] += 1

        inventory.append({
            "dataset": dataset.upper(),
            "partition": filename,
            **dict(counts),
            "unique_direct_hla_pairs": len(file_pairs)
        })

        print(
            f"{filename}: {counts['rows']:,} rows; "
            f"{counts['direct_hla_rows']:,} direct HLA rows; "
            f"{counts['cell_line_rows']:,} cell-line rows"
        )

    # Compare reference data with IEDB.

    exact = iedb_pairs & direct_pairs
    possible = iedb_pairs & possible_pairs
    sequence_only = iedb_sequences & reference_sequences

    for peptide, allele in sorted(exact):
        all_exact_overlaps.append({
            "reference_dataset": dataset.upper(),
            "peptide": peptide,
            "hla_allele": allele,
            "overlap_type": "Direct allele assignment"
        })

    for peptide, allele in sorted(possible):
        if (peptide, allele) not in exact:
            all_possible_overlaps.append({
                "reference_dataset": dataset.upper(),
                "peptide": peptide,
                "hla_allele": allele,
                "overlap_type": "Multi-allelic possible association"
            })

    overlap_summary.append({
        "reference_dataset": dataset.upper(),
        "reference_direct_pairs": len(direct_pairs),
        "reference_possible_MA_pairs": len(possible_pairs),
        "direct_exact_pair_overlaps": len(exact),
        "possible_MA_pair_overlaps": len(possible - exact),
        "sequence_only_overlaps": len(sequence_only)
    })

    print("Direct exact overlaps:", len(exact))
    print("Possible multi-allelic overlaps:", len(possible - exact))
    print("Sequence-only overlaps:", len(sequence_only))


# Save results.


def save_csv(filename, rows):

    if not rows:
        print(f"Skipped {filename}: no records to save.")
        return

    # Collect every unique column across all rows.
    fieldnames = list(dict.fromkeys(
        key
        for row in rows
        for key in row.keys()
    ))

    with (OUTPUT / filename).open(
        "w", newline="", encoding="utf-8-sig"
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            restval=""
        )

        writer.writeheader()

        for row in rows:
            clean_row = {}

            for key, value in row.items():

                # Preserve structured values as valid JSON strings.
                if isinstance(value, (dict, list, tuple)):
                    value = json.dumps(
                        value,
                        ensure_ascii=False
                    )

                clean_row[key] = value

            writer.writerow(clean_row)

    print(
        f"Saved {filename}: "
        f"{len(rows):,} rows, {len(fieldnames)} columns"
    )

print("\n========== AUDIT COMPLETE ==========")
print("Inventory:", OUTPUT / "netmhcpan_training_inventory.csv")
print("Overlap summary:", OUTPUT / "netmhcpan_iedb_overlap_pre_split.csv")
print("Exact pairs:", OUTPUT / "netmhcpan_exact_overlap_pairs.csv")
print("Possible MA pairs:", OUTPUT / "netmhcpan_possible_MA_overlap_pairs.csv")
