"""End-to-end wiring tests for hseq-importer lifecycle.

Validates the complete chain from MainService through HseqImportService
to PostgresHseqRepository, ensuring the lifecycle order is correct,
run_comment/run_submitted propagate consistently, and Genew4Lock is
used for pointer updates.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch, call

import pytest

from hgnc_hseq_importer.models import HseqCandidate
from hgnc_hseq_importer.services.hseq_import_service import HseqImportService
from hgnc_hseq_importer.services.main_service import MainService


class TestEndToEndLifecycle:
    """Tests for the full MainService → HseqImportService → Repository chain."""

    def test_full_lifecycle_order(self) -> None:
        mock_repo = MagicMock()
        mock_repo.get_pseudogene_candidates.return_value = []
        mock_repo.get_vega_candidates.return_value = []
        mock_repo.get_ccds_candidates.return_value = []
        mock_repo.get_ensembl_candidates.return_value = []

        lock = MagicMock()
        service = HseqImportService(
            repository=mock_repo,
            logger=MagicMock(),
            genew4_lock=lock,
        )

        service.run_import()

        assert mock_repo.method_calls[0].args == ()
        assert mock_repo.batch_insert_hseq.call_count == 0
        assert mock_repo.update_hgnc_hseq_pointers.call_count == 0

    def test_lifecycle_with_candidates_propagates_params(self) -> None:
        candidates = [
            HseqCandidate(
                hgnc_id=1,
                source="ccds",
                defline="CCDS1 | C:1 | EG:100 | SYM | HGNC:1",
                sequence="ATCG",
                status="bulk",
            ),
        ]

        mock_repo = MagicMock()
        mock_repo.get_pseudogene_candidates.return_value = []
        mock_repo.get_vega_candidates.return_value = []
        mock_repo.get_ccds_candidates.return_value = candidates
        mock_repo.get_ensembl_candidates.return_value = []
        mock_repo.batch_insert_hseq.return_value = 1
        mock_repo.update_hgnc_hseq_pointers.return_value = 1

        lock = MagicMock()
        service = HseqImportService(
            repository=mock_repo,
            logger=MagicMock(),
            genew4_lock=lock,
        )

        result = service.run_import()

        insert_call = mock_repo.batch_insert_hseq.call_args
        assert insert_call.kwargs["run_comment"] == "import via hseqs_importer"
        insert_submitted = insert_call.kwargs["run_submitted"]

        pointer_call = mock_repo.update_hgnc_hseq_pointers.call_args
        assert pointer_call.kwargs["run_comment"] == "import via hseqs_importer"
        assert pointer_call.kwargs["run_submitted"] == insert_submitted
        assert pointer_call.kwargs["genew4_lock"] is lock

        assert result.sequences_upserted == 1
        assert result.pointers_updated == 1

    def test_from_settings_wires_genew4_lock_to_service(self) -> None:
        mock_lock_module = MagicMock()
        mock_lock_instance = MagicMock()
        mock_lock_module.Genew4Lock.return_value = mock_lock_instance
        mock_lock_module.Genew4LockSqlRepository.return_value = MagicMock()

        mock_session = MagicMock()

        with patch("hgnc_hseq_importer.services.main_service.create_engine"):
            with patch(
                "hgnc_hseq_importer.services.main_service.Session",
                return_value=mock_session,
            ):
                with patch.dict("sys.modules", {"shared.genew4_lock": mock_lock_module}):
                    with patch.dict(
                        "sys.modules",
                        {"shared": MagicMock(genew4_lock=mock_lock_module)},
                    ):
                        settings = MagicMock()
                        settings.genew4.dsn.return_value = (
                            "postgresql://user:pass@host:5432/db"
                        )
                        service = MainService.from_settings(settings)

        assert service._import_service._genew4_lock is mock_lock_instance
