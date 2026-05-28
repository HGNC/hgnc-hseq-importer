"""Abstract repository interface for hseq and HGNC pointer operations.

Defines the contract that concrete Postgres repository implementations
must fulfil for reading/writing sequence records and updating HGNC
sequence pointer columns.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from hgnc_hseq_importer.models import HgncGene, HseqRecord


class HseqRepository(ABC):
    """Define the persistence contract for hseq sequence storage.

    All methods are batch-oriented. Implementations must ensure
    idempotency: re-running with the same input yields no net changes.
    """

    @abstractmethod
    def get_sequences_for_gene(self, hgnc_id: str) -> list[HseqRecord]:
        """Return all hseq records for the given gene.

        Args:
            hgnc_id: The HGNC identifier.

        Returns:
            A list of HseqRecord instances.
        """

    @abstractmethod
    def upsert_sequences(self, records: list[HseqRecord]) -> int:
        """Insert or update hseq records.

        Args:
            records: The sequence records to upsert.

        Returns:
            The number of rows affected.
        """

    @abstractmethod
    def get_hgnc_genes_with_sequences(self) -> list[HgncGene]:
        """Return all HGNC genes that have at least one sequence.

        Returns:
            A list of HgncGene instances with current pointer values.
        """

    @abstractmethod
    def update_hgnc_pointers(self, genes: list[HgncGene]) -> int:
        """Update HGNC sequence pointer columns for the given genes.

        Args:
            genes: Genes with updated pointer values to persist.

        Returns:
            The number of rows updated.
        """
