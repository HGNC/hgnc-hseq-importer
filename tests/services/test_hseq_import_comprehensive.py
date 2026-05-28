"""Comprehensive tests for priority chain and end-to-end flow (Task 36.8).

Verify the full coordinator flow from source queries through priority
selection, including edge cases from the Perl reference behaviour.
"""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from hgnc_hseq_importer.models import HseqCandidate
from hgnc_hseq_importer.services.hseq_import_service import (
    HseqImportCoordinator,
    HseqImportService,
)


class TestFullPriorityChain:
    """Test the full priority chain: pseudo → vega → ccds → ensembl."""

    def test_pseudo_gene_skips_vega_ccds_ensembl(self) -> None:
        coordinator = HseqImportCoordinator()

        pseudo = HseqCandidate(
            hgnc_id=1, source="pseudo",
            defline="BRCA1 | 12345 | C:17 | HGNC:1", sequence="A"
        )
        vega = HseqCandidate(
            hgnc_id=1, source="vega",
            defline="OTTHUMG001 | HGNC:1", sequence="B"
        )
        ccds = HseqCandidate(
            hgnc_id=1, source="ccds",
            defline="CCDS1 | C:17 | EG:672 | BRCA1 | HGNC:1", sequence="C"
        )
        ensembl = HseqCandidate(
            hgnc_id=1, source="ensembl",
            defline="ENSG000001 | HGNC:1", sequence="D"
        )

        result = coordinator.select_candidates([pseudo, vega, ccds, ensembl])

        assert len(result) == 1
        assert result[0].source == "pseudo"

    def test_vega_gene_skips_ccds_and_ensembl(self) -> None:
        coordinator = HseqImportCoordinator()

        vega = HseqCandidate(
            hgnc_id=5, source="vega",
            defline="OTTHUMG005 | HGNC:5", sequence="B"
        )
        ccds = HseqCandidate(
            hgnc_id=5, source="ccds",
            defline="CCDS5 | C:1 | EG:100 | SYM | HGNC:5", sequence="C"
        )
        ensembl = HseqCandidate(
            hgnc_id=5, source="ensembl",
            defline="ENSG000005 | HGNC:5", sequence="D"
        )

        result = coordinator.select_candidates([vega, ccds, ensembl])

        assert len(result) == 1
        assert result[0].source == "vega"

    def test_ccds_gene_skips_ensembl(self) -> None:
        coordinator = HseqImportCoordinator()

        ccds = HseqCandidate(
            hgnc_id=10, source="ccds",
            defline="CCDS10 | C:2 | EG:200 | SYM | HGNC:10", sequence="C"
        )
        ensembl = HseqCandidate(
            hgnc_id=10, source="ensembl",
            defline="ENSG000010 | HGNC:10", sequence="D"
        )

        result = coordinator.select_candidates([ccds, ensembl])

        assert len(result) == 1
        assert result[0].source == "ccds"

    def test_ensembl_gene_when_no_higher_priority(self) -> None:
        coordinator = HseqImportCoordinator()

        ensembl = HseqCandidate(
            hgnc_id=20, source="ensembl",
            defline="ENSG000020 | HGNC:20", sequence="D"
        )

        result = coordinator.select_candidates([ensembl])

        assert len(result) == 1
        assert result[0].source == "ensembl"


class TestMultipleGenesWithMixedSources:
    """Test realistic scenarios with multiple genes from different sources."""

    def test_mixed_sources_across_genes(self) -> None:
        coordinator = HseqImportCoordinator()

        candidates = [
            HseqCandidate(hgnc_id=1, source="pseudo", defline="d", sequence="A"),
            HseqCandidate(hgnc_id=2, source="vega", defline="d", sequence="B"),
            HseqCandidate(hgnc_id=3, source="ccds", defline="d", sequence="C"),
            HseqCandidate(hgnc_id=4, source="ensembl", defline="d", sequence="D"),
        ]

        result = coordinator.select_candidates(candidates)

        assert len(result) == 4
        sources = {c.hgnc_id: c.source for c in result}
        assert sources[1] == "pseudo"
        assert sources[2] == "vega"
        assert sources[3] == "ccds"
        assert sources[4] == "ensembl"

    def test_gene_without_any_source_not_in_result(self) -> None:
        coordinator = HseqImportCoordinator()

        candidates = [
            HseqCandidate(hgnc_id=1, source="pseudo", defline="d", sequence="A"),
        ]

        result = coordinator.select_candidates(candidates)

        assert all(c.hgnc_id != 99 for c in result)

    def test_large_batch_priority_correctness(self) -> None:
        coordinator = HseqImportCoordinator()

        candidates = []
        for i in range(100):
            candidates.append(
                HseqCandidate(
                    hgnc_id=i, source="pseudo", defline=f"d{i}", sequence=f"seq{i}"
                )
            )
            candidates.append(
                HseqCandidate(
                    hgnc_id=i, source="ensembl", defline=f"d{i}", sequence=f"seq{i}"
                )
            )

        result = coordinator.select_candidates(candidates)

        assert len(result) == 100
        assert all(c.source == "pseudo" for c in result)


class TestCoordinatorStateTracking:
    """Test that coordinator correctly tracks handled genes across calls."""

    def test_handled_ids_persist_across_select_calls(self) -> None:
        coordinator = HseqImportCoordinator()

        coordinator.select_candidates([
            HseqCandidate(hgnc_id=1, source="pseudo", defline="d", sequence="A"),
        ])

        result = coordinator.select_candidates([
            HseqCandidate(hgnc_id=1, source="vega", defline="d", sequence="B"),
            HseqCandidate(hgnc_id=2, source="vega", defline="d", sequence="C"),
        ])

        assert len(result) == 1
        assert result[0].hgnc_id == 2

    def test_handled_hgnc_ids_returns_all(self) -> None:
        coordinator = HseqImportCoordinator()

        coordinator.select_candidates([
            HseqCandidate(hgnc_id=10, source="pseudo", defline="d", sequence="A"),
            HseqCandidate(hgnc_id=20, source="vega", defline="d", sequence="B"),
        ])

        assert coordinator.handled_hgnc_ids() == {10, 20}


class TestHseqImportServiceOrchestration:
    """Test HseqImportService.run_import with mocked repository."""

    def test_run_import_queries_all_sources(self) -> None:
        mock_repo = MagicMock()
        mock_repo.get_pseudogene_candidates.return_value = []
        mock_repo.get_vega_candidates.return_value = []
        mock_repo.get_ccds_candidates.return_value = []
        mock_repo.get_ensembl_candidates.return_value = []
        mock_repo.batch_insert_hseq.return_value = 0

        import logging

        service = HseqImportService(
            repository=mock_repo,
            logger=logging.getLogger("test"),
        )
        result = service.run_import()

        mock_repo.get_pseudogene_candidates.assert_called_once()
        mock_repo.get_vega_candidates.assert_called_once()
        mock_repo.get_ccds_candidates.assert_called_once()
        mock_repo.get_ensembl_candidates.assert_called_once()

    def test_run_import_inserts_priority_selected_candidates(self) -> None:
        mock_repo = MagicMock()
        mock_repo.get_pseudogene_candidates.return_value = [
            HseqCandidate(hgnc_id=1, source="pseudo", defline="d1", sequence="A"),
        ]
        mock_repo.get_vega_candidates.return_value = [
            HseqCandidate(hgnc_id=1, source="vega", defline="d2", sequence="B"),
            HseqCandidate(hgnc_id=2, source="vega", defline="d3", sequence="C"),
        ]
        mock_repo.get_ccds_candidates.return_value = []
        mock_repo.get_ensembl_candidates.return_value = []
        mock_repo.batch_insert_hseq.return_value = 2

        import logging

        service = HseqImportService(
            repository=mock_repo,
            logger=logging.getLogger("test"),
        )
        result = service.run_import()

        insert_args = mock_repo.batch_insert_hseq.call_args[0][0]
        inserted_ids = {c.hgnc_id for c in insert_args}
        assert 1 in inserted_ids
        assert 2 in inserted_ids
        assert len(insert_args) == 2
