"""Sequence parser for multi-source sequence ingestion.

Parses raw TSV content into validated HseqRecord domain models,
normalizes identifiers, and handles malformed rows gracefully.
"""

from __future__ import annotations

import logging

from hgnc_hseq_importer.exceptions import ParseError
from hgnc_hseq_importer.models import HseqRecord

logger = logging.getLogger(__name__)


class SequenceParser:
    """Parse raw TSV sequence input into validated HseqRecord models.

    Expects tab-separated fields:
    hgnc_id, source, sequence_type, sequence, accession[, version]

    Malformed rows are logged and skipped.
    """

    def parse(self, raw: str) -> list[HseqRecord]:
        """Parse raw TSV content into a list of HseqRecord instances.

        Args:
            raw: The raw TSV content.

        Returns:
            A list of validated HseqRecord instances.

        Raises:
            ParseError: If the input is empty.
        """
        if not raw or not raw.strip():
            raise ParseError("Sequence input is empty")

        records: list[HseqRecord] = []
        skipped = 0

        for line_no, line in enumerate(raw.strip().splitlines(), start=1):
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                logger.warning("Skipping malformed row %d: insufficient columns", line_no)
                skipped += 1
                continue

            try:
                version = int(parts[5]) if len(parts) > 5 else 1
                record = HseqRecord(
                    hgnc_id=parts[0].strip(),
                    source=parts[1].strip(),
                    sequence_type=parts[2].strip(),
                    sequence=parts[3].strip(),
                    accession=parts[4].strip(),
                    version=version,
                )
                records.append(record)
            except (ValueError, Exception) as exc:
                logger.warning("Skipping invalid row %d: %s", line_no, exc)
                skipped += 1

        return records


def normalize_record(record: HseqRecord) -> HseqRecord:
    """Normalize a single HseqRecord by stripping whitespace.

    Args:
        record: The raw record.

    Returns:
        A normalized record with trimmed fields.
    """
    return HseqRecord(
        hgnc_id=record.hgnc_id.strip(),
        source=record.source.strip(),
        sequence_type=record.sequence_type.strip(),
        sequence=record.sequence.strip(),
        accession=record.accession.strip(),
        version=record.version,
    )


def normalize_batch(records: list[HseqRecord]) -> list[HseqRecord]:
    """Normalize a batch of HseqRecord instances.

    Args:
        records: The raw records.

    Returns:
        A list of normalized records.
    """
    return [normalize_record(r) for r in records]
