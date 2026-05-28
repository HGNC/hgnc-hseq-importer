"""Domain models for hseq sequence storage and HGNC pointer updates.

Defines Pydantic models representing sequence records, HGNC gene
pointers, and a source precedence enum for deterministic canonical
sequence selection.
"""

from __future__ import annotations

import enum

from pydantic import BaseModel, Field


class SourcePrecedence(str, enum.Enum):
    """Rank sources for canonical sequence selection.

    Lower value = higher precedence. When multiple sequences exist for
    the same gene and sequence_type, the source with the lowest rank is
    chosen as canonical. Ties within the same source are broken by
    accession lexicographic order for determinism.
    """

    ENSEMBL = "1"
    REFSEQ = "2"
    UNIPROT = "3"


class HseqRecord(BaseModel):
    """Represent a single sequence record from an external source.

    Attributes:
        hgnc_id: The HGNC identifier linking this sequence to a gene.
        source: The source database (e.g. "ensembl", "refseq").
        sequence_type: The type of sequence (e.g. "peptide", "cdna").
        sequence: The raw sequence string.
        accession: The source-specific accession identifier.
        version: Optional version number for the sequence.
    """

    hgnc_id: str = Field(min_length=1)
    source: str = Field(min_length=1)
    sequence_type: str = Field(min_length=1)
    sequence: str = Field(default="")
    accession: str = Field(min_length=1)
    version: int = Field(default=1)


class HgncGene(BaseModel):
    """Represent an HGNC gene with its current sequence pointer columns.

    Attributes:
        hgnc_id: The HGNC identifier.
        symbol: The approved gene symbol.
        peptide_refseq_id: Current RefSeq peptide accession pointer (or None).
        peptide_ensembl_id: Current Ensembl peptide accession pointer (or None).
    """

    hgnc_id: str = Field(min_length=1)
    symbol: str = Field(default="")
    peptide_refseq_id: str | None = None
    peptide_ensembl_id: str | None = None
