"""Tests for hseq-importer MainService wiring.

Validates that MainService correctly constructs sessions,
PostgresHseqRepository, Genew4Lock, and HseqImportService from settings,
and that run() executes the full import lifecycle.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from hgnc_hseq_importer.services.main_service import MainService


def _make_settings() -> MagicMock:
    s = MagicMock()
    s.genew4.dsn.return_value = "postgresql://user:pass@host:5432/db"
    return s


class TestMainServiceWiring:
    """Tests for MainService dependency wiring."""

    @patch("hgnc_hseq_importer.services.main_service.Session")
    @patch("hgnc_hseq_importer.services.main_service.create_engine")
    def test_from_settings_creates_engine(
        self, mock_engine_factory, mock_session_cls
    ) -> None:
        mock_lock_module = MagicMock()
        with patch(
            "hgnc_hseq_importer.services.main_service.genew4_lock",
            mock_lock_module,
            create=True,
        ):
            with patch.dict("sys.modules", {"shared.genew4_lock": mock_lock_module}):
                with patch.dict("sys.modules", {"shared": MagicMock(genew4_lock=mock_lock_module)}):
                    MainService.from_settings(_make_settings())

        mock_engine_factory.assert_called_once()

    @patch("hgnc_hseq_importer.services.main_service.Session")
    @patch("hgnc_hseq_importer.services.main_service.create_engine")
    def test_from_settings_creates_two_sessions(
        self, mock_engine_factory, mock_session_cls
    ) -> None:
        mock_lock_module = MagicMock()
        with patch.dict("sys.modules", {"shared.genew4_lock": mock_lock_module}):
            with patch.dict("sys.modules", {"shared": MagicMock(genew4_lock=mock_lock_module)}):
                MainService.from_settings(_make_settings())

        assert mock_session_cls.call_count >= 2

    @patch("hgnc_hseq_importer.services.main_service.Session")
    @patch("hgnc_hseq_importer.services.main_service.create_engine")
    def test_from_settings_creates_lock_with_name(
        self, mock_engine_factory, mock_session_cls
    ) -> None:
        mock_lock_module = MagicMock()
        with patch.dict("sys.modules", {"shared.genew4_lock": mock_lock_module}):
            with patch.dict("sys.modules", {"shared": MagicMock(genew4_lock=mock_lock_module)}):
                MainService.from_settings(_make_settings())

        mock_lock_module.Genew4Lock.assert_called_once()
        call_kwargs = mock_lock_module.Genew4Lock.call_args
        assert call_kwargs.kwargs.get("name") == "hseqs_importer" or call_kwargs[1].get("name") == "hseqs_importer"

    def test_run_delegates_to_import_service(self) -> None:
        mock_import_service = MagicMock()
        service = MainService(import_service=mock_import_service)
        service.run()

        mock_import_service.run.assert_called_once()
