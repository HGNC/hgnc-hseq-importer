"""Domain models for hseq sequence storage and HGNC pointer updates.

Defines the HseqCandidate dataclass used by the coordinator and
repository layers. The old SourcePrecedence enum, HseqRecord,
and HgncGene models have been removed in favour of the
coordinator-based priority algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class HseqCandidate:
    """A candidate Hseq record from a source query.

    Attributes:
        hgnc_id: The HGNC identifier for the gene.
        source: The source name (pseudo, vega, ccds, ensembl).
        defline: The FASTA defline string.
        sequence: The nucleotide sequence.
    """

    hgnc_id: int
    source: str
    defline: str
    sequence: str
