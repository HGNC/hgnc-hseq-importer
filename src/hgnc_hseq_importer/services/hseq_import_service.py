"""HseqImportService orchestrating sequence ingestion and pointer recomputation.

Coordinates the full import flow: parses raw input, normalizes sequences,
upserts them via the repository, selects canonical sequences, recomputes
HGNC pointers, and applies the updates.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass

from hgnc_hseq_importer.canonical_selector import (
    recompute_pointers,
    select_canonical_sequences,
)
from hgnc_hseq_importer.models import HgncGene
from hgnc_hseq_importer.repositories.hseq_repository import HseqRepository
from hgnc_hseq_importer.sequence_parser import normalize_batch

if __name__ == "__main__":  # pragma: no cover
    pass
else:
    from hgnc_hseq_importer.sequence_parser import SequenceParser


@dataclass
class ImportResult:
    """Structured summary of an hseq import run.

    Attributes:
        records_parsed: Total records parsed from input.
        sequences_upserted: Sequences inserted or updated.
        pointers_updated: HGNC pointer columns changed.
        pointers_unchanged: HGNC pointer columns unchanged.
    """

    records_parsed: int = 0
    sequences_upserted: int = 0
    pointers_updated: int = 0
    pointers_unchanged: int = 0


class HseqImportService:
    """Orchestrate the hseq import and HGNC pointer recomputation lifecycle.

    Args:
        parser: Sequence parser for converting raw input to records.
        repository: Repository for persisting sequences and updating pointers.
        logger: Logger for emitting structured metrics.
    """

    def __init__(
        self,
        parser: SequenceParser,
        repository: HseqRepository,
        logger: logging.Logger,
    ) -> None:
        self._parser = parser
        self._repository = repository
        self._logger = logger

    def run_import(self, raw_input: str) -> ImportResult:
        """Execute the full hseq import lifecycle.

        Args:
            raw_input: The raw TSV input data.

        Returns:
            An ImportResult with counts of parsed, upserted, and updated records.

        Raises:
            PersistenceError: If database operations fail.
        """
        self._logger.info(
            "import_start",
            extra={"event": "import_start"},
        )

        start = time.monotonic()
        try:
            parsed = self._parser.parse(raw_input)
            normalized = normalize_batch(parsed)
            upserted = self._repository.upsert_sequences(normalized)

            genes = self._repository.get_hgnc_genes_with_sequences()
            canonical = select_canonical_sequences(normalized)
            updated_genes = recompute_pointers(canonical, genes)

            genes_with_changes = [
                g for g in updated_genes
                if self._pointer_changed(g, next(
                    (og for og in genes if og.hgnc_id == g.hgnc_id), None
                ))
            ]
            pointers_updated = 0
            if genes_with_changes:
                pointers_updated = self._repository.update_hgnc_pointers(genes_with_changes)

            result = ImportResult(
                records_parsed=len(parsed),
                sequences_upserted=upserted,
                pointers_updated=pointers_updated,
                pointers_unchanged=len(updated_genes) - pointers_updated,
            )
        finally:
            elapsed = time.monotonic() - start

        self._logger.info(
            "import_complete",
            extra={
                "event": "import_complete",
                "records_parsed": result.records_parsed,
                "sequences_upserted": result.sequences_upserted,
                "pointers_updated": result.pointers_updated,
                "pointers_unchanged": result.pointers_unchanged,
                "duration_seconds": round(elapsed, 3),
            },
        )

        return result

    def _pointer_changed(self, new: HgncGene, old: HgncGene | None) -> bool:
        """Check if a gene's pointers have changed.

        Args:
            new: The recomputed gene with updated pointers.
            old: The original gene, or None if new.

        Returns:
            True if any pointer column differs.
        """
        if old is None:
            return True
        return (
            new.peptide_ensembl_id != old.peptide_ensembl_id
            or new.peptide_refseq_id != old.peptide_refseq_id
        )
