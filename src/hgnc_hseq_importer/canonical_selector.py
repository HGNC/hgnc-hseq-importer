"""Canonical sequence selection and HGNC pointer recomputation.

Implements deterministic canonical selection using source precedence and
accession lexicographic tie-breaking. Recomputes HGNC sequence pointer
columns based on the selected canonical sequences.
"""

from __future__ import annotations

from collections import defaultdict

from hgnc_hseq_importer.models import HgncGene, HseqRecord, SourcePrecedence


def _source_rank(source: str) -> int:
    """Return the precedence rank for a source name.

    Lower is better. Unknown sources get the lowest precedence.

    Args:
        source: The source name.

    Returns:
        An integer rank for sorting.
    """
    source_upper = source.upper()
    for member in SourcePrecedence:
        if member.name == source_upper:
            return int(member.value)
    return 999


def select_canonical_sequences(sequences: list[HseqRecord]) -> list[HseqRecord]:
    """Select one canonical peptide sequence per gene using source precedence.

    For each gene, picks the sequence with the lowest source precedence rank.
    Ties within the same source are broken by accession lexicographic order
    for determinism. Only considers sequences of type "peptide".

    Args:
        sequences: All available sequence records.

    Returns:
        One canonical sequence per gene that has peptide records.
    """
    if not sequences:
        return []

    by_gene: dict[str, list[HseqRecord]] = defaultdict(list)
    for seq in sequences:
        if seq.sequence_type == "peptide":
            by_gene[seq.hgnc_id].append(seq)

    canonical: list[HseqRecord] = []
    for hgnc_id, gene_seqs in by_gene.items():
        gene_seqs.sort(key=lambda s: (_source_rank(s.source), s.accession))
        canonical.append(gene_seqs[0])

    return canonical


def recompute_pointers(
    canonical: list[HseqRecord],
    current_genes: list[HgncGene],
) -> list[HgncGene]:
    """Recompute HGNC sequence pointer columns from canonical sequences.

    For each gene, sets the peptide_ensembl_id and peptide_refseq_id
    pointers based on the canonical sequences. Pointers not backed by
    a canonical sequence are set to None.

    Args:
        canonical: The selected canonical sequences.
        current_genes: Current HGNC gene records with existing pointers.

    Returns:
        Updated HgncGene instances with recomputed pointer values.
    """
    canonical_by_gene: dict[str, dict[str, str]] = defaultdict(dict)
    for seq in canonical:
        canonical_by_gene[seq.hgnc_id][seq.source.lower()] = seq.accession

    updated: list[HgncGene] = []
    for gene in current_genes:
        pointers = canonical_by_gene.get(gene.hgnc_id, {})
        updated.append(
            HgncGene(
                hgnc_id=gene.hgnc_id,
                symbol=gene.symbol,
                peptide_ensembl_id=pointers.get("ensembl"),
                peptide_refseq_id=pointers.get("refseq"),
            )
        )

    return updated
