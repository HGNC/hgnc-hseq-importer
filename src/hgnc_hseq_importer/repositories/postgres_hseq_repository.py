"""Postgres repository for hseq and HGNC pointer operations.

Implements the HseqRepository ABC using psycopg v3 with parameterized
SQL for all operations.
"""

from __future__ import annotations

import logging
from typing import Any

from hgnc_hseq_importer.exceptions import PersistenceError
from hgnc_hseq_importer.models import HgncGene, HseqRecord
from hgnc_hseq_importer.repositories.hseq_repository import HseqRepository

logger = logging.getLogger(__name__)

_SELECT_SEQUENCES_FOR_GENE = (
    "SELECT hgnc_id, source, sequence_type, sequence, accession, version "
    "FROM hseq WHERE hgnc_id = %s"
)

_UPSERT_SEQUENCE = (
    "INSERT INTO hseq (hgnc_id, source, sequence_type, sequence, accession, version) "
    "VALUES (%s, %s, %s, %s, %s, %s) "
    "ON CONFLICT (hgnc_id, source, sequence_type, accession) DO UPDATE SET "
    "sequence = EXCLUDED.sequence, "
    "version = EXCLUDED.version"
)

_SELECT_GENES_WITH_SEQUENCES = (
    "SELECT g.hgnc_id, g.symbol, "
    "g.peptide_ensembl_id, g.peptide_refseq_id "
    "FROM hgnc g WHERE EXISTS ("
    "  SELECT 1 FROM hseq s WHERE s.hgnc_id = g.hgnc_id"
    ")"
)

_UPDATE_HGNC_POINTERS = (
    "UPDATE hgnc SET "
    "peptide_ensembl_id = %s, "
    "peptide_refseq_id = %s "
    "WHERE hgnc_id = %s"
)


class PostgresHseqRepository(HseqRepository):
    """Concrete Postgres repository for hseq sequence storage.

    Uses psycopg v3 with parameterized SQL. Maps database errors
    to PersistenceError for service-layer consumption.

    Args:
        connection: A psycopg v3 connection instance.
    """

    def __init__(self, connection: Any) -> None:
        self._conn = connection

    def get_sequences_for_gene(self, hgnc_id: str) -> list[HseqRecord]:
        try:
            with self._conn.cursor() as cur:
                cur.execute(_SELECT_SEQUENCES_FOR_GENE, (hgnc_id,))
                return [
                    HseqRecord(
                        hgnc_id=row[0],
                        source=row[1],
                        sequence_type=row[2],
                        sequence=row[3],
                        accession=row[4],
                        version=row[5],
                    )
                    for row in cur.fetchall()
                ]
        except Exception as exc:
            raise PersistenceError(
                f"Failed to query hseq for gene {hgnc_id!r}: {exc}"
            ) from exc

    def upsert_sequences(self, records: list[HseqRecord]) -> int:
        if not records:
            return 0
        try:
            params = [
                (r.hgnc_id, r.source, r.sequence_type, r.sequence, r.accession, r.version)
                for r in records
            ]
            with self._conn.cursor() as cur:
                cur.executemany(_UPSERT_SEQUENCE, params)
                return cur.rowcount
        except Exception as exc:
            raise PersistenceError(
                f"Failed to upsert {len(records)} hseq records: {exc}"
            ) from exc

    def get_hgnc_genes_with_sequences(self) -> list[HgncGene]:
        try:
            with self._conn.cursor() as cur:
                cur.execute(_SELECT_GENES_WITH_SEQUENCES)
                return [
                    HgncGene(
                        hgnc_id=row[0],
                        symbol=row[1],
                        peptide_ensembl_id=row[2],
                        peptide_refseq_id=row[3],
                    )
                    for row in cur.fetchall()
                ]
        except Exception as exc:
            raise PersistenceError(
                f"Failed to query HGNC genes with sequences: {exc}"
            ) from exc

    def update_hgnc_pointers(self, genes: list[HgncGene]) -> int:
        if not genes:
            return 0
        try:
            params = [
                (g.peptide_ensembl_id, g.peptide_refseq_id, g.hgnc_id)
                for g in genes
            ]
            with self._conn.cursor() as cur:
                cur.executemany(_UPDATE_HGNC_POINTERS, params)
                return cur.rowcount
        except Exception as exc:
            raise PersistenceError(
                f"Failed to update HGNC pointers for {len(genes)} genes: {exc}"
            ) from exc
