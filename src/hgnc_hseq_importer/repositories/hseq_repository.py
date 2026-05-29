"""Abstract repository interface for hseq and HGNC pointer operations.

Defines the contract that concrete repository implementations must
fulfil for querying source data and persisting Hseq records.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from hgnc_hseq_importer.models import HseqCandidate


class HseqRepository(ABC):
    """Define the persistence contract for hseq sequence storage.

    All methods are batch-oriented. Implementations must ensure
    idempotency: re-running with the same input yields no net changes.
    """

    @abstractmethod
    def get_pseudogene_candidates(self) -> list[HseqCandidate]:
        """Return Pseudogene source candidates.

        Returns:
            A list of HseqCandidate instances from pseudogene_org.
        """

    @abstractmethod
    def get_vega_candidates(self) -> list[HseqCandidate]:
        """Return VEGA source candidates.

        Returns:
            A list of HseqCandidate instances from otter_seq.
        """

    @abstractmethod
    def get_ccds_candidates(self) -> list[HseqCandidate]:
        """Return CCDS source candidates.

        Returns:
            A list of HseqCandidate instances from ccds + ccds_seq.
        """

    @abstractmethod
    def get_ensembl_candidates(self) -> list[HseqCandidate]:
        """Return Ensembl source candidates.

        Returns:
            A list of HseqCandidate instances from ensembl_seq.
        """

    @abstractmethod
    def batch_insert_hseq(self, candidates: list[HseqCandidate]) -> int:
        """Insert Hseq records for the given candidates.

        Args:
            candidates: The candidate records to persist.

        Returns:
            The number of rows inserted.
        """

    @abstractmethod
    def update_hgnc_hseq_pointers(
        self,
        run_comment: str,
        run_submitted: int,
        editor: str,
        genew4_lock: object | None = None,
    ) -> int:
        """Update Gene.hseq_ids and Gene.pub_hseq_id for newly inserted Hseq rows.

        Finds new hseq rows matching the run metadata, parses hgnc_id
        from the defline, and updates the corresponding Gene records
        using Genew4Lock for safe concurrent access.

        Args:
            run_comment: The comment used to identify this run's hseq rows.
            run_submitted: The submitted timestamp for this run.
            editor: The editor name for lock acquisition.
            genew4_lock: A Genew4Lock instance for row-level locking.

        Returns:
            The number of Gene records updated.
        """
