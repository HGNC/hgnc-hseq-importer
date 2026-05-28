"""Main service for the HGNC HSEQ importer."""

from typing import TYPE_CHECKING

from hgnc_hseq_importer.services.base_service import Service

if TYPE_CHECKING:
    from hgnc_hseq_importer.config import Settings


class MainService(Service):
    """Orchestrate the HGNC HSEQ importing workflow.

    This service coordinates the steps required to import HSEQ
    data. All database access is delegated to injected repositories.
    """

    def __init__(self, settings: "Settings") -> None:
        self._settings = settings

    @classmethod
    def from_settings(cls, settings: "Settings") -> "MainService":
        """Construct a MainService from a Settings instance.

        Args:
            settings: Application configuration.

        Returns:
            A configured MainService ready to run.
        """
        return cls(settings=settings)

    def run(self) -> None:
        """Execute the HSEQ importing workflow."""
