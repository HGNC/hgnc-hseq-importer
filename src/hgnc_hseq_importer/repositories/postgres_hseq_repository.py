"""Postgres repository for hseq and HGNC pointer operations.

Implements the HseqRepository ABC using genew4-orm SQLAlchemy sessions
for source queries and psycopg v3 for the Hseq insert / Gene update
lifecycle.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from genew4_orm import models as orm
from sqlalchemy import select, func

from hgnc_hseq_importer.models import HseqCandidate
from hgnc_hseq_importer.exceptions import PersistenceError
from hgnc_hseq_importer.repositories.hseq_repository import HseqRepository

if TYPE_CHECKING:
    from sqlalchemy.orm import Session


class PostgresHseqRepository(HseqRepository):
    """Concrete Postgres repository for hseq sequence storage.

    Uses genew4-orm SQLAlchemy sessions for queries. Maps database
    errors to PersistenceError for service-layer consumption.

    Args:
        readonly_session: A SQLAlchemy Session for read queries.
        readwrite_session: A SQLAlchemy Session for writes.
    """

    def __init__(
        self,
        readonly_session: Session,
        readwrite_session: Session,
    ) -> None:
        self._ro = readonly_session
        self._rw = readwrite_session

    def get_pseudogene_candidates(self) -> list[HseqCandidate]:
        try:
            stmt = (
                select(
                    orm.Gene.hgnc_id,
                    orm.Gene.approved_symbol,
                    orm.Gene.pseudogene_id,
                    orm.PseudogeneOrg.sequence,
                    orm.PseudogeneOrg.chromosome,
                )
                .select_from(orm.Gene, orm.PseudogeneOrg)
                .where(orm.Gene.pseudogene_id == orm.PseudogeneOrg.porg_id)
                .where(orm.Gene.hseq_ids.is_(None))
                .where(orm.Gene.lock.is_(None))
            )
            rows = self._ro.execute(stmt).all()
            return [
                HseqCandidate(
                    hgnc_id=row[0],
                    source="pseudo",
                    defline=f"{row[1]} | {row[2]} | C:{row[4]} | HGNC:{row[0]}",
                    sequence=row[3] or "",
                )
                for row in rows
            ]
        except Exception as exc:
            raise PersistenceError(
                f"Failed to query pseudogene candidates: {exc}"
            ) from exc

    def get_vega_candidates(self) -> list[HseqCandidate]:
        try:
            stmt = (
                select(
                    orm.Gene.hgnc_id,
                    orm.OtterSequence.oseq_gene_id,
                    orm.OtterSequence.defline,
                    orm.OtterSequence.sequence,
                )
                .select_from(orm.Gene, orm.OtterSequence)
                .where(orm.Gene.vega_ids == orm.OtterSequence.oseq_gene_id)
                .where(orm.Gene.pseudogene_id.is_(None))
                .where(orm.Gene.hseq_ids.is_(None))
                .where(orm.Gene.lock.is_(None))
                .where(orm.Gene.public_refseq_ids.is_(None))
                .order_by(orm.OtterSequence.length.desc())
            )
            rows = self._ro.execute(stmt).all()
            seen_gene_ids: set[str] = set()
            candidates: list[HseqCandidate] = []
            for row in rows:
                gene_id = row[1]
                if gene_id in seen_gene_ids:
                    continue
                seen_gene_ids.add(gene_id)
                candidates.append(
                    HseqCandidate(
                        hgnc_id=row[0],
                        source="vega",
                        defline=f"{row[2]} | HGNC:{row[0]}",
                        sequence=row[3] or "",
                    )
                )
            return candidates
        except Exception as exc:
            raise PersistenceError(
                f"Failed to query VEGA candidates: {exc}"
            ) from exc

    def get_ccds_candidates(self) -> list[HseqCandidate]:
        try:
            first_ccds = func.split_part(orm.Gene.ccds_ids, ",", 1)
            stmt = (
                select(
                    orm.Gene.hgnc_id,
                    orm.Ccds.ccds_id,
                    orm.Ccds.chromosome,
                    orm.Ccds.ncbi_gene_id,
                    orm.Gene.approved_symbol,
                    orm.CcdsSequence.sequence,
                )
                .select_from(orm.Gene)
                .join(orm.Ccds, first_ccds == orm.Ccds.ccds_id)
                .join(orm.CcdsSequence, orm.Ccds.ccds_id == orm.CcdsSequence.ccdseq_ccds_id)
                .where(orm.Gene.hseq_ids.is_(None))
                .where(orm.Gene.public_refseq_ids.is_(None))
            )
            rows = self._ro.execute(stmt).all()
            return [
                HseqCandidate(
                    hgnc_id=row[0],
                    source="ccds",
                    defline=(
                        f"{row[1]} | C:{row[2]} | EG:{row[3]} | {row[4]} | HGNC:{row[0]}"
                    ),
                    sequence=row[5] or "",
                    status="bulk",
                )
                for row in rows
            ]
        except Exception as exc:
            raise PersistenceError(
                f"Failed to query CCDS candidates: {exc}"
            ) from exc

    def get_ensembl_candidates(self) -> list[HseqCandidate]:
        try:
            stmt = (
                select(
                    orm.Gene.hgnc_id,
                    orm.EnsemblSequence.defline,
                    orm.EnsemblSequence.sequence,
                )
                .select_from(orm.Gene, orm.EnsemblSequence)
                .where(orm.Gene.public_ensembl_id == orm.EnsemblSequence.eseq_ensembl_gene_id)
                .where(orm.Gene.hseq_ids.is_(None))
                .where(orm.Gene.public_refseq_ids.is_(None))
                .order_by(
                    orm.EnsemblSequence.length.desc(),
                    orm.EnsemblSequence.eseq_ensembl_gene_id,
                )
            )
            rows = self._ro.execute(stmt).all()
            return [
                HseqCandidate(
                    hgnc_id=row[0],
                    source="ensembl",
                    defline=f"{row[1]} | HGNC:{row[0]}",
                    sequence=row[2] or "",
                )
                for row in rows
            ]
        except Exception as exc:
            raise PersistenceError(
                f"Failed to query Ensembl candidates: {exc}"
            ) from exc

    def batch_insert_hseq(self, candidates: list[HseqCandidate]) -> int:
        if not candidates:
            return 0
        try:
            records = [
                orm.Hseq(
                    ext=c.source,
                    editor="genew",
                    molecule="dna",
                    submitted=int(__import__("time").time()),
                    status=c.status,
                    priority=100,
                    run_notes="search=hgnc_heavy, summ=50, align=30",
                    comment="import via hseqs_importer",
                    entry_class="archive",
                    is_new="TRUE",
                    defline=c.defline,
                    sequence=c.sequence,
                )
                for c in candidates
            ]
            self._rw.add_all(records)
            self._rw.flush()
            return len(records)
        except Exception as exc:
            raise PersistenceError(
                f"Failed to batch insert {len(candidates)} hseq records: {exc}"
            ) from exc

    def update_hgnc_hseq_pointers(
        self,
        run_comment: str,
        run_submitted: int,
        editor: str,
    ) -> int:
        raise NotImplementedError("TODO: implement with Genew4Lock")
