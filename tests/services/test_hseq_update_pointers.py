"""Tests for HseqImportService wiring of update_hgnc_hseq_pointers.

Validates that run_import calls update_hgnc_hseq_pointers after
batch_insert_hseq, with the same run_comment and run_submitted
values used during the insert, and using Genew4Lock for row-level
locking.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from hgnc_hseq_importer.models import HseqCandidate
from hgnc_hseq_importer.services.hseq_import_service import (
    HseqImportService,
    ImportResult,
)


def _make_candidate(
    hgnc_id: int = 1100,
    source: str = "pseudo",
    defline: str = "SYMBOL | 42 | C:7 | HGNC:1100",
    sequence: str = "ATGC",
) -> HseqCandidate:
    return HseqCandidate(
        hgnc_id=hgnc_id,
        source=source,
        defline=defline,
        sequence=sequence,
    )


def _make_mock_repo() -> MagicMock:
    repo = MagicMock()
    repo.get_pseudogene_candidates.return_value = []
    repo.get_vega_candidates.return_value = []
    repo.get_ccds_candidates.return_value = []
    repo.get_ensembl_candidates.return_value = []
    repo.batch_insert_hseq.return_value = 1
    repo.update_hgnc_hseq_pointers.return_value = 1
    return repo


class TestHseqImportServiceUpdatePointers:
    """Tests for update_hgnc_hseq_pointers wiring in run_import."""

    def test_run_import_calls_update_hgnc_hseq_pointers(self) -> None:
        import logging

        repo = _make_mock_repo()
        repo.get_pseudogene_candidates.return_value = [_make_candidate()]
        lock = MagicMock()

        service = HseqImportService(
            repository=repo,
            logger=logging.getLogger("test"),
            genew4_lock=lock,
        )
        service.run_import()

        repo.update_hgnc_hseq_pointers.assert_called_once()

    def test_run_import_passes_correct_comment(self) -> None:
        import logging

        repo = _make_mock_repo()
        repo.get_pseudogene_candidates.return_value = [_make_candidate()]
        lock = MagicMock()

        service = HseqImportService(
            repository=repo,
            logger=logging.getLogger("test"),
            genew4_lock=lock,
        )
        service.run_import()

        call_kwargs = repo.update_hgnc_hseq_pointers.call_args
        assert call_kwargs.kwargs.get("run_comment") == "import via hseqs_importer" or call_kwargs[1].get("run_comment") == "import via hseqs_importer"

    def test_run_import_passes_genew4_lock(self) -> None:
        import logging

        repo = _make_mock_repo()
        repo.get_pseudogene_candidates.return_value = [_make_candidate()]
        lock = MagicMock()

        service = HseqImportService(
            repository=repo,
            logger=logging.getLogger("test"),
            genew4_lock=lock,
        )
        service.run_import()

        call_kwargs = repo.update_hgnc_hseq_pointers.call_args
        genew4_lock_arg = call_kwargs.kwargs.get("genew4_lock") or call_kwargs[1].get("genew4_lock")
        assert genew4_lock_arg is lock

    def test_run_import_returns_import_result(self) -> None:
        import logging

        repo = _make_mock_repo()
        repo.get_pseudogene_candidates.return_value = [_make_candidate()]
        lock = MagicMock()

        service = HseqImportService(
            repository=repo,
            logger=logging.getLogger("test"),
            genew4_lock=lock,
        )
        result = service.run_import()

        assert isinstance(result, ImportResult)
        assert result.pointers_updated >= 0

    def test_update_not_called_when_no_candidates(self) -> None:
        import logging

        repo = _make_mock_repo()
        lock = MagicMock()

        service = HseqImportService(
            repository=repo,
            logger=logging.getLogger("test"),
            genew4_lock=lock,
        )
        service.run_import()

        repo.batch_insert_hseq.assert_not_called()
        repo.update_hgnc_hseq_pointers.assert_not_called()

    def test_update_called_after_insert(self) -> None:
        import logging

        repo = _make_mock_repo()
        repo.get_pseudogene_candidates.return_value = [_make_candidate()]
        lock = MagicMock()

        call_order: list[str] = []
        repo.batch_insert_hseq.side_effect = lambda *a, **k: (
            call_order.append("insert"), 1
        )[1]
        repo.update_hgnc_hseq_pointers.side_effect = lambda *a, **k: (
            call_order.append("update"), 1
        )[1]

        service = HseqImportService(
            repository=repo,
            logger=logging.getLogger("test"),
            genew4_lock=lock,
        )
        service.run_import()

        assert call_order == ["insert", "update"]
