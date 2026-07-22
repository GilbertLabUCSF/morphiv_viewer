#!/usr/bin/env python3
"""Build a filterable Parquet copy of the EBs pathway-enrichment export."""

from __future__ import annotations

from pathlib import Path

import pyarrow.csv as csv
import pyarrow.parquet as parquet


PROJECT_ROOT = Path(__file__).resolve().parent.parent
SOURCE = (
    PROJECT_ROOT
    / "data"
    / "tf_perturbseq"
    / "results"
    / "EBs"
    / "pathway_enrichment"
    / "EBs_enrichment_results.csv"
)
OUTPUT = PROJECT_ROOT / "data_extracted" / "pathway_enrichment_ebs.parquet"


def main() -> int:
    if not SOURCE.exists():
        print(f"Pathway enrichment source not found: {SOURCE}")
        return 1

    OUTPUT.parent.mkdir(exist_ok=True)
    temporary = OUTPUT.with_suffix(".parquet.tmp")
    reader = csv.open_csv(
        SOURCE,
        read_options=csv.ReadOptions(block_size=32 * 1024 * 1024),
    )

    writer = None
    row_count = 0
    try:
        for batch in reader:
            if writer is None:
                writer = parquet.ParquetWriter(
                    temporary,
                    batch.schema,
                    compression="zstd",
                    use_dictionary=True,
                    write_statistics=True,
                )
            # Small row groups let the portal use perturbation min/max statistics
            # to skip almost the entire table for a single-gene query.
            writer.write_batch(batch, row_group_size=2_000)
            row_count += batch.num_rows
    finally:
        if writer is not None:
            writer.close()

    if row_count == 0:
        temporary.unlink(missing_ok=True)
        print("Pathway enrichment source is empty")
        return 1

    temporary.replace(OUTPUT)
    print(
        f"Wrote {row_count:,} pathway rows to {OUTPUT} "
        f"({OUTPUT.stat().st_size / 1e6:.1f} MB)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
