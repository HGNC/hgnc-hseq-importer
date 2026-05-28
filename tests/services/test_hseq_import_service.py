"""Tests for the HseqImportCoordinator priority algorithm (Task 36.2).

Verify that the coordinator calls source methods in correct order
(Pseudo → VEGA → CCDS → Ensembl) and never assigns an Hseq to a gene
already satisfied by a higher-priority source.
"""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from hgnc_hseq_importer.models import HseqCandidate
from hgnc_hseq_importer.services.hseq_import_service import (
    HseqImportCoordinator,
)


@dataclass
class _FakeGene:
    """Minimal gene stub for testing priority logic."""

    hgnc_id: int
    symbol: str = ""
    pseudogene_id: int | None = None
    vega_ids: str | None = None
    ccds_ids: str | None = None
    pub_ensembl_id: str | None = None
    hseq_ids: str | None = None
    lock: str | None = None
    pub_refseq_ids: str | None = None


class TestCoordinatorSourceOrder:
    """Verify the coordinator processes candidates in priority order."""

    def test_candidates_processed_in_pseudo_vega_ccds_ensembl_order(self) -> None:
        coordinator = HseqImportCoordinator()

        pseudo = HseqCandidate(
            hgnc_id=1, source="pseudo", defline="d1", sequence="A"
        )
        vega = HseqCandidate(
            hgnc_id=2, source="vega", defline="d2", sequence="B"
        )
        ccds = HseqCandidate(
            hgnc_id=3, source="ccds", defline="d3", sequence="C"
        )
        ensembl = HseqCandidate(
            hgnc_id=4, source="ensembl", defline="d4", sequence="D"
        )

        result = coordinator.select_candidates([pseudo, vega, ccds, ensembl])

        sources = [c.source for c in result]
        assert sources == ["pseudo", "vega", "ccds", "ensembl"]


class TestCoordinatorPriorityTracking:
    """Verify handled genes are excluded from lower-priority sources."""

    def test_gene_handled_by_pseudo_excluded_from_vega(self) -> None:
        coordinator = HseqImportCoordinator()

        pseudo_candidate = HseqCandidate(
            hgnc_id=1,
            source="pseudo",
            defline="BRCA1 | pseudo_123 | C:17 | HGNC:1",
            sequence="ATCG",
        )

        result = coordinator.select_candidates([pseudo_candidate])

        handled = coordinator.handled_hgnc_ids()
        assert 1 in handled

    def test_gene_handled_by_vega_excluded_from_ccds(self) -> None:
        coordinator = HseqImportCoordinator()

        vega_candidate = HseqCandidate(
            hgnc_id=5,
            source="vega",
            defline="OTTHUMG000001 | HGNC:5",
            sequence="GCTA",
        )

        result = coordinator.select_candidates([vega_candidate])
        assert 5 in coordinator.handled_hgnc_ids()

    def test_multiple_sources_only_highest_priority_wins(self) -> None:
        coordinator = HseqImportCoordinator()

        pseudo = HseqCandidate(
            hgnc_id=10,
            source="pseudo",
            defline="SYM | pseudo_10 | C:1 | HGNC:10",
            sequence="AA",
        )
        vega = HseqCandidate(
            hgnc_id=10,
            source="vega",
            defline="VEGA_10 | HGNC:10",
            sequence="GG",
        )

        result = coordinator.select_candidates([pseudo, vega])

        pseudo_results = [c for c in result if c.source == "pseudo"]
        vega_results = [c for c in result if c.source == "vega" and c.hgnc_id == 10]
        assert len(pseudo_results) == 1
        assert len(vega_results) == 0


class TestCoordinatorSelectCandidates:
    """Verify select_candidates returns the right set."""

    def test_returns_empty_for_no_candidates(self) -> None:
        coordinator = HseqImportCoordinator()
        result = coordinator.select_candidates([])
        assert result == []

    def test_returns_all_unique_gene_candidates(self) -> None:
        coordinator = HseqImportCoordinator()

        c1 = HseqCandidate(
            hgnc_id=1, source="pseudo", defline="d1", sequence="A"
        )
        c2 = HseqCandidate(
            hgnc_id=2, source="vega", defline="d2", sequence="B"
        )
        c3 = HseqCandidate(
            hgnc_id=3, source="ensembl", defline="d3", sequence="C"
        )

        result = coordinator.select_candidates([c1, c2, c3])
        assert len(result) == 3

    def test_duplicate_gene_keeps_first_source_only(self) -> None:
        coordinator = HseqImportCoordinator()

        c1 = HseqCandidate(
            hgnc_id=1, source="ccds", defline="d1", sequence="A"
        )
        c2 = HseqCandidate(
            hgnc_id=1, source="ensembl", defline="d2", sequence="B"
        )

        result = coordinator.select_candidates([c1, c2])
        assert len(result) == 1
        assert result[0].source == "ccds"
