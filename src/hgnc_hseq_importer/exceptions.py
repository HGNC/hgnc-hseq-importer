"""Domain-specific exception hierarchy for the HGNC HSEQ importer."""


class ServiceError(Exception):
    """Base exception for all domain errors."""


class ConfigError(ServiceError):
    """Raised when configuration is invalid or missing."""


class RepositoryError(ServiceError):
    """Raised when a database or data access operation fails."""


class ParseError(ServiceError):
    """Raised when parsing input data fails."""


class PersistenceError(ServiceError):
    """Raised when a database write or update operation fails."""


class ValidationError(ServiceError):
    """Raised when domain validation of a record fails."""
