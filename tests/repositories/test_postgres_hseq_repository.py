"""Tests for Pseudogene and VEGA source methods (Task 36.4-36.5).

Verify source query methods return correct candidates with proper join
conditions, filters, deduplication, and defline formatting.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from hgnc_hseq_importer.models import HseqCandidate
from hgnc_hseq_importer.repositories.postgres_hseq_repository import (
    PostgresHseqRepository,
)


class TestGetPseudogeneCandidates:
    """Test get_pseudogene_candidates query and defline formatting."""

    def test_returns_candidates_with_correct_defline(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = [
            (1, "BRCA1", 12345, "ATCGATCG", "17"),
            (2, "TP53", 67890, "GCTAGCTA", "17"),
        ]
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_pseudogene_candidates()

        assert len(candidates) == 2
        assert candidates[0].hgnc_id == 1
        assert candidates[0].source == "pseudo"
        assert candidates[0].defline == "BRCA1 | 12345 | C:17 | HGNC:1"
        assert candidates[0].sequence == "ATCGATCG"
        assert candidates[1].defline == "TP53 | 67890 | C:17 | HGNC:2"

    def test_handles_null_sequence(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = [
            (5, "SYM", 999, None, "1"),
        ]
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_pseudogene_candidates()

        assert candidates[0].sequence == ""

    def test_returns_empty_for_no_matches(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = []
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_pseudogene_candidates()

        assert candidates == []

    def test_raises_persistence_error_on_db_failure(self) -> None:
        from hgnc_hseq_importer.exceptions import PersistenceError

        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_ro.execute.side_effect = RuntimeError("connection lost")

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )

        with pytest.raises(PersistenceError, match="pseudogene candidates"):
            repo.get_pseudogene_candidates()


class TestGetVegaCandidates:
    """Test get_vega_candidates with deduplication and defline formatting."""

    def test_returns_candidates_with_correct_defline(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = [
            (10, "OTTHUMG000001", "OTTHUMG000001-001 some desc", "AAAACCCC", None),
        ]
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_vega_candidates()

        assert len(candidates) == 1
        assert candidates[0].hgnc_id == 10
        assert candidates[0].source == "vega"
        assert candidates[0].defline == "OTTHUMG000001-001 some desc | HGNC:10"
        assert candidates[0].sequence == "AAAACCCC"

    def test_deduplicates_by_gene_id(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = [
            (20, "OTTHUMG000002", "longest desc", "AAAACCCCGGGG", None),
            (20, "OTTHUMG000002", "shorter desc", "AA", None),
        ]
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_vega_candidates()

        assert len(candidates) == 1
        assert candidates[0].defline == "longest desc | HGNC:20"

    def test_handles_null_sequence(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = [
            (30, "OTTHUMG000003", "desc", None, None),
        ]
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_vega_candidates()

        assert candidates[0].sequence == ""

    def test_returns_empty_for_no_matches(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = []
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_vega_candidates()

        assert candidates == []

    def test_raises_persistence_error_on_db_failure(self) -> None:
        from hgnc_hseq_importer.exceptions import PersistenceError

        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_ro.execute.side_effect = RuntimeError("connection lost")

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )

        with pytest.raises(PersistenceError, match="VEGA candidates"):
            repo.get_vega_candidates()


class TestGetCcdsCandidates:
    """Test get_ccds_candidates with first-CCDS and defline formatting."""

    def test_returns_candidates_with_correct_defline(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = [
            (50, "CCDS1", "17", "1234", "BRCA1", "ATCGATCG"),
        ]
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_ccds_candidates()

        assert len(candidates) == 1
        assert candidates[0].hgnc_id == 50
        assert candidates[0].source == "ccds"
        assert candidates[0].defline == "CCDS1 | C:17 | EG:1234 | BRCA1 | HGNC:50"
        assert candidates[0].sequence == "ATCGATCG"

    def test_handles_null_sequence(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = [
            (51, "CCDS2", "1", "5678", "SYM", None),
        ]
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_ccds_candidates()

        assert candidates[0].sequence == ""

    def test_returns_empty_for_no_matches(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = []
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_ccds_candidates()

        assert candidates == []

    def test_raises_persistence_error_on_db_failure(self) -> None:
        from hgnc_hseq_importer.exceptions import PersistenceError

        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_ro.execute.side_effect = RuntimeError("connection lost")

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )

        with pytest.raises(PersistenceError, match="CCDS candidates"):
            repo.get_ccds_candidates()


class TestGetEnsemblCandidates:
    """Test get_ensembl_candidates with defline formatting."""

    def test_returns_candidates_with_correct_defline(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = [
            (100, "ENSG000001 desc", "GGGGCCCC"),
        ]
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_ensembl_candidates()

        assert len(candidates) == 1
        assert candidates[0].hgnc_id == 100
        assert candidates[0].source == "ensembl"
        assert candidates[0].defline == "ENSG000001 desc | HGNC:100"
        assert candidates[0].sequence == "GGGGCCCC"

    def test_handles_null_sequence(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = [
            (101, "some desc", None),
        ]
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_ensembl_candidates()

        assert candidates[0].sequence == ""

    def test_returns_empty_for_no_matches(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_result = MagicMock()
        mock_result.all.return_value = []
        mock_ro.execute.return_value = mock_result

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        candidates = repo.get_ensembl_candidates()

        assert candidates == []

    def test_raises_persistence_error_on_db_failure(self) -> None:
        from hgnc_hseq_importer.exceptions import PersistenceError

        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_ro.execute.side_effect = RuntimeError("connection lost")

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )

        with pytest.raises(PersistenceError, match="Ensembl candidates"):
            repo.get_ensembl_candidates()


class TestHseqCandidateStatus:
    """Test HseqCandidate status field defaults and overrides."""

    def test_default_status_is_done(self) -> None:
        candidate = HseqCandidate(
            hgnc_id=1, source="pseudo", defline="d", sequence="A"
        )
        assert candidate.status == "done"

    def test_status_can_be_overridden_to_bulk(self) -> None:
        candidate = HseqCandidate(
            hgnc_id=1, source="ccds", defline="d", sequence="A", status="bulk"
        )
        assert candidate.status == "bulk"

    def test_pseudo_source_uses_default_status(self) -> None:
        candidate = HseqCandidate(
            hgnc_id=1, source="pseudo", defline="d", sequence="A"
        )
        assert candidate.status == "done"

    def test_vega_source_uses_default_status(self) -> None:
        candidate = HseqCandidate(
            hgnc_id=1, source="vega", defline="d", sequence="A"
        )
        assert candidate.status == "done"

    def test_ensembl_source_uses_default_status(self) -> None:
        candidate = HseqCandidate(
            hgnc_id=1, source="ensembl", defline="d", sequence="A"
        )
        assert candidate.status == "done"


class TestBatchInsertHseq:
    """Test batch_insert_hseq with candidate records."""

    def test_inserts_records_and_returns_count(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()

        candidates = [
            HseqCandidate(
                hgnc_id=1, source="pseudo", defline="d1", sequence="A"
            ),
            HseqCandidate(
                hgnc_id=2, source="vega", defline="d2", sequence="B"
            ),
        ]

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        count = repo.batch_insert_hseq(candidates)

        assert count == 2
        mock_rw.add_all.assert_called_once()
        mock_rw.flush.assert_called_once()

    def test_returns_zero_for_empty_list(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        count = repo.batch_insert_hseq([])

        assert count == 0
        mock_rw.add_all.assert_not_called()

    def test_raises_persistence_error_on_db_failure(self) -> None:
        from hgnc_hseq_importer.exceptions import PersistenceError

        mock_ro = MagicMock()
        mock_rw = MagicMock()
        mock_rw.add_all.side_effect = RuntimeError("disk full")

        candidates = [
            HseqCandidate(
                hgnc_id=1, source="pseudo", defline="d1", sequence="A"
            ),
        ]

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )

        with pytest.raises(PersistenceError, match="batch insert"):
            repo.batch_insert_hseq(candidates)

    def test_ccds_candidates_use_bulk_status(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        added_records = []
        mock_rw.add_all.side_effect = lambda records: added_records.extend(records)

        candidates = [
            HseqCandidate(
                hgnc_id=1, source="ccds", defline="d1", sequence="A",
                status="bulk",
            ),
            HseqCandidate(
                hgnc_id=2, source="pseudo", defline="d2", sequence="B",
            ),
            HseqCandidate(
                hgnc_id=3, source="vega", defline="d3", sequence="C",
            ),
            HseqCandidate(
                hgnc_id=4, source="ensembl", defline="d4", sequence="D",
            ),
        ]

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        repo.batch_insert_hseq(candidates)

        assert len(added_records) == 4
        assert added_records[0].status == "bulk"
        assert added_records[1].status == "done"
        assert added_records[2].status == "done"
        assert added_records[3].status == "done"

    def test_non_ccds_candidates_use_done_status(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()
        added_records = []
        mock_rw.add_all.side_effect = lambda records: added_records.extend(records)

        candidates = [
            HseqCandidate(
                hgnc_id=1, source="pseudo", defline="d1", sequence="A",
            ),
            HseqCandidate(
                hgnc_id=2, source="vega", defline="d2", sequence="B",
            ),
            HseqCandidate(
                hgnc_id=3, source="ensembl", defline="d3", sequence="C",
            ),
        ]

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )
        repo.batch_insert_hseq(candidates)

        assert all(r.status == "done" for r in added_records)


class TestUpdateHgncHseqPointers:
    """Test update_hgnc_hseq_pointers (currently NotImplementedError)."""

    def test_raises_not_implemented(self) -> None:
        mock_ro = MagicMock()
        mock_rw = MagicMock()

        repo = PostgresHseqRepository(
            readonly_session=mock_ro, readwrite_session=mock_rw
        )

        with pytest.raises(NotImplementedError):
            repo.update_hgnc_hseq_pointers(
                run_comment="test",
                run_submitted=12345,
                editor="genew",
            )
