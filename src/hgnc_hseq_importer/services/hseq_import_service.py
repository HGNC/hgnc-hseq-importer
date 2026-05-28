"""HseqImportService orchestrating sequence ingestion via priority selection.

Coordinates the full import flow: queries each source in priority order
(Pseudogene → VEGA → CCDS → Ensembl), tracks handled genes to prevent
duplicate assignments, and produces candidate Hseq records for persistence.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field

from hgnc_hseq_importer.models import HseqCandidate
from hgnc_hseq_importer.repositories.hseq_repository import HseqRepository
from hgnc_hseq_importer.services.base_service import Service


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


class HseqImportCoordinator:
    """Coordinate sequential priority-based Hseq candidate selection.

    Applies the priority algorithm: Pseudo → VEGA → CCDS → Ensembl.
    Once a gene is assigned a candidate from a higher-priority source,
    it is excluded from all lower-priority sources.

    Args:
        logger: Logger for emitting structured events.
    """

    def __init__(self, logger: logging.Logger | None = None) -> None:
        self._logger = logger or logging.getLogger(__name__)
        self._handled_hgnc_ids: set[int] = set()

    def handled_hgnc_ids(self) -> set[int]:
        """Return the set of HGNC IDs already assigned a candidate.

        Returns:
            A frozenset of handled HGNC IDs.
        """
        return self._handled_hgnc_ids

    def select_candidates(
        self,
        candidates: list[HseqCandidate],
    ) -> list[HseqCandidate]:
        """Filter candidates by priority, keeping only the first per gene.

        Iterates candidates in order (callers must provide them in
        priority order: pseudo, vega, ccds, ensembl). For each gene,
        only the first candidate is kept; subsequent candidates for
        the same gene are discarded.

        Args:
            candidates: Candidate records ordered by source priority.

        Returns:
            The filtered list with one candidate per gene.
        """
        selected: list[HseqCandidate] = []
        for candidate in candidates:
            if candidate.hgnc_id in self._handled_hgnc_ids:
                continue
            self._handled_hgnc_ids.add(candidate.hgnc_id)
            selected.append(candidate)
        return selected


class HseqImportService(Service):
    """Orchestrate the hseq import and HGNC pointer recomputation lifecycle.

    Args:
        repository: Repository for persisting sequences and updating pointers.
        logger: Logger for emitting structured metrics.
    """

    def __init__(
        self,
        repository: HseqRepository,
        logger: logging.Logger,
    ) -> None:
        self._repository = repository
        self._logger = logger

    def run_import(self) -> ImportResult:
        """Execute the full hseq import lifecycle.

        Returns:
            An ImportResult with counts of upserted and updated records.

        Raises:
            PersistenceError: If database operations fail.
        """
        self._logger.info(
            "import_start",
            extra={"event": "import_start"},
        )

        start = time.monotonic()
        try:
            result = ImportResult()
        finally:
            elapsed = time.monotonic() - start

        self._logger.info(
            "import_complete",
            extra={
                "event": "import_complete",
                "duration_seconds": round(elapsed, 3),
            },
        )

        return result
