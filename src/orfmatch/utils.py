from __future__ import annotations

from typing import Dict, List, Tuple
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def log(message: str) -> None:
    print(f"[orfmatch] {message}")


def parse_protein_id(record: SeqRecord) -> str:
    """Return a stable locus/id from a FASTA record header.

    Uses the first whitespace-delimited token of the FASTA header (record.id).
    """
    return str(record.id).split()[0]


def parse_protein_fasta(
    reference_faa: str,
) -> Tuple[List[SeqRecord], List[SeqRecord], Dict[str, SeqFeature], Dict[str, List[SeqFeature]]]:
    """Parse a protein FASTA and build reference lists and feature maps.

    Returns: (reference_proteins, reference_rnas, protein_feature_map, rna_feature_map)
      - reference_rnas and rna_feature_map are empty for FASTA input.
    """
    reference_proteins: List[SeqRecord] = []
    reference_rnas: List[SeqRecord] = []  # no RNAs in protein FASTA
    protein_feature_map: Dict[str, SeqFeature] = {}
    rna_feature_map: Dict[str, List[SeqFeature]] = {}

    for record in SeqIO.parse(reference_faa, "fasta"):
        locus = parse_protein_id(record)
        seq_str = str(record.seq)
        protein = SeqRecord(Seq(seq_str), id=locus, description="")
        reference_proteins.append(protein)

        # Minimal feature so downstream code (Annotator) can copy qualifiers
        feat = SeqFeature(
            FeatureLocation(0, len(seq_str), strand=1),
            type="CDS",
            qualifiers={
                "translation": [seq_str],
                "locus_tag": [locus],
                "product": [str(record.id)],
            },
        )
        protein_feature_map[locus] = feat

    return reference_proteins, reference_rnas, protein_feature_map, rna_feature_map


def parse_gbff(
    reference_gbff: str,
) -> Tuple[List[SeqRecord], List[SeqRecord], Dict[str, SeqFeature], Dict[str, List[SeqFeature]]]:
    """Parse a GenBank/GBFF file and build reference lists and feature maps.

    Returns: (reference_proteins, reference_rnas, protein_feature_map, rna_feature_map)
    """
    reference_proteins: List[SeqRecord] = []
    reference_rnas: List[SeqRecord] = []
    protein_feature_map: Dict[str, SeqFeature] = {}
    rna_feature_map: Dict[str, List[SeqFeature]] = {}

    rna_id_counter = defaultdict(int)

    for record in SeqIO.parse(reference_gbff, "genbank"):
        for feature in getattr(record, "features", []):
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                prot_seq = str(feature.qualifiers["translation"][0])
                prot_seq = prot_seq.rstrip("*")
                locus = feature.qualifiers.get("locus_tag", [record.id])[0]
                protein = SeqRecord(Seq(prot_seq), id=locus, description="")
                reference_proteins.append(protein)
                protein_feature_map[locus] = feature
            elif feature.type in {"tRNA", "rRNA", "ncRNA"}:
                rna_seq = feature.extract(record.seq)
                locus = feature.qualifiers.get("locus_tag", [record.id])[0]
                rna_id_counter[locus] += 1
                # ensure uniqueness if multiple RNAs share the same locus_tag
                unique_id = locus if rna_id_counter[
                    locus] == 1 else f"{locus}_{rna_id_counter[locus]}"
                reference_rnas.append(
                    SeqRecord(rna_seq, id=unique_id, description=feature.type))
                rna_feature_map.setdefault(unique_id, []).append(feature)

    return reference_proteins, reference_rnas, protein_feature_map, rna_feature_map


__all__ = [
    "log",
    "parse_protein_id",
    "parse_protein_fasta",
    "parse_gbff",
]
