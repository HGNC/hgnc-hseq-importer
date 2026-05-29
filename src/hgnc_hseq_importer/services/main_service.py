"""Main orchestrator for the HGNC HSEQ importer.

Constructs all dependencies from Settings, wires the full lifecycle,
and delegates to HseqImportService for execution.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from hgnc_hseq_importer.repositories.postgres_hseq_repository import (
    PostgresHseqRepository,
)
from hgnc_hseq_importer.services.base_service import Service
from hgnc_hseq_importer.services.hseq_import_service import HseqImportService

if TYPE_CHECKING:
    from hgnc_hseq_importer.config import Settings

logger = logging.getLogger(__name__)


class MainService(Service):
    """Orchestrate the HGNC HSEQ importing workflow.

    Wires up the SQLAlchemy engine, sessions, Genew4Lock, repository,
    and HseqImportService. All database access is delegated to the
    injected repository.

    Args:
        import_service: The HseqImportService that performs the actual
            import lifecycle.
    """

    def __init__(self, import_service: HseqImportService) -> None:
        """Initialise with a fully wired HseqImportService.

        Args:
            import_service: Service that handles the import lifecycle.
        """
        self._import_service = import_service

    @classmethod
    def from_settings(cls, settings: Settings) -> MainService:
        """Construct a MainService from a Settings instance.

        Creates the SQLAlchemy engine with pool resilience, opens
        separate read-only and read-write sessions, builds the
        Genew4Lock, and wires everything into HseqImportService.

        Args:
            settings: Application configuration with genew4 DSN.

        Returns:
            A configured MainService ready to run.
        """
        from shared import genew4_lock

        engine = create_engine(
            settings.genew4.dsn(),
            pool_pre_ping=True,
            pool_recycle=300,
        )

        ro_session = Session(engine)
        rw_session = Session(engine)

        repository = PostgresHseqRepository(
            readonly_session=ro_session,
            readwrite_session=rw_session,
        )

        lock_repo = genew4_lock.Genew4LockSqlRepository(session=rw_session)
        lock = genew4_lock.Genew4Lock(
            repository=lock_repo,
            name="hseqs_importer",
        )

        import_service = HseqImportService(
            repository=repository,
            logger=logger,
            genew4_lock=lock,
        )

        return cls(import_service=import_service)

    def run(self) -> None:
        """Execute the HSEQ importing workflow.

        Delegates to HseqImportService.run_import() and logs the result.
        """
        self._import_service.run()
