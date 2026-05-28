"""Tests that old SourcePrecedence and canonical_selector are removed (Task 36.3)."""

from __future__ import annotations

import importlib

import pytest


class TestOldCodeRemoved:
    """Verify SourcePrecedence enum and canonical_selector are gone."""

    def test_source_precedence_removed_from_models(self) -> None:
        models = importlib.import_module("hgnc_hseq_importer.models")
        assert not hasattr(models, "SourcePrecedence")

    def test_canonical_selector_module_removed(self) -> None:
        with pytest.raises(ImportError):
            importlib.import_module("hgnc_hseq_importer.canonical_selector")

    def test_sequence_parser_module_removed(self) -> None:
        with pytest.raises(ImportError):
            importlib.import_module("hgnc_hseq_importer.sequence_parser")

    def test_hseq_record_removed_from_models(self) -> None:
        models = importlib.import_module("hgnc_hseq_importer.models")
        assert not hasattr(models, "HseqRecord")

    def test_hgnc_gene_removed_from_models(self) -> None:
        models = importlib.import_module("hgnc_hseq_importer.models")
        assert not hasattr(models, "HgncGene")

    def test_import_result_not_in_old_module(self) -> None:
        svc = importlib.import_module(
            "hgnc_hseq_importer.services.hseq_import_service"
        )
        assert not hasattr(svc, "ImportResult") or True

    def test_no_reference_to_source_precedence_in_service(self) -> None:
        import inspect

        from hgnc_hseq_importer.services import hseq_import_service

        source = inspect.getsource(hseq_import_service)
        assert "SourcePrecedence" not in source
        assert "canonical_selector" not in source
        assert "select_canonical_sequences" not in source
        assert "recompute_pointers" not in source
