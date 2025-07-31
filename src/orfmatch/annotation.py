from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

import pyrodigal
from pyhmmer import easel, hmmer
from tqdm import tqdm

from orfmatch.utils import log

import os


class Annotator:
    """
    A callable annotation pipeline you can import and use from a separate
    `main.py`. Construct it with your inputs, then call `.annotate()`.

    Example usage from main.py:
        from orfmatch.annotation import Annotator
        annot = Annotator(
            contigs, reference_proteins, reference_rnas,
            reference_gbff="refs.gbff",
            protein_feature_map=protein_feature_map,
            rna_feature_map=rna_feature_map,
            annotated_gbff="annotated.gbff",
            variants_fasta="variants.faa",
            threads=8, e_value=1e-5, show_variants=True,
        )
        records = annot.annotate()
    """

    def __init__(
        self,
        contigs: List[SeqRecord],
        reference_proteins: Iterable[SeqRecord],
        reference_rnas: Iterable[SeqRecord],
        reference_gbff: str,
        *,
        protein_feature_map: Dict[str, SeqFeature],
        rna_feature_map: Dict[str, List[SeqFeature]],
        annotated_gbff: str = "annotated.gbff",
        variants_fasta: str = "variants.faa",
        threads: int = 4,
        e_value: float = 1e-5,
        show_variants: bool = False,
    ) -> None:
        # Inputs
        self.contigs = contigs
        self.reference_proteins = list(reference_proteins)
        self.reference_rnas = list(reference_rnas)
        self.reference_gbff = reference_gbff
        self.protein_feature_map = protein_feature_map
        self.rna_feature_map = rna_feature_map

        # Settings
        self.annotated_gbff = annotated_gbff
        self.variants_fasta = variants_fasta
        self.threads = threads
        self.e_value = e_value
        self.show_variants = show_variants

        # Alphabets
        self.prot_alphabet = easel.Alphabet.amino()
        self.dna_alphabet = easel.Alphabet.dna()

        # Precompute digital references for HMMER and an exact-match lookup
        self.digital_refs: Dict[str, easel.DigitalSequence] = {
            rec.id: easel.TextSequence(name=rec.id.encode(), sequence=str(rec.seq)).digitize(self.prot_alphabet)
            for rec in self.reference_proteins
        }
        self.exact_ref_lookup: Dict[str, str] = {
            str(rec.seq).rstrip("*").strip(): rec.id for rec in self.reference_proteins
        }
        self.digital_rnas: Dict[str, List[easel.DigitalSequence]] = defaultdict(list)
        for rec in self.reference_rnas:
            digital = easel.TextSequence(name=rec.id.encode(), sequence=str(rec.seq)).digitize(self.dna_alphabet)
            self.digital_rnas[rec.id].append(digital)

    # -------------------- ORF finding --------------------
    def _find_orfs(self) -> Tuple[List[SeqRecord], List[Tuple[SeqFeature, str]]]:
        """Predict ORFs on the contigs using pyrodigal.
        Returns (prodigal_records, predicted_features) where predicted_features is
        a list of (CDS SeqFeature, translation str).
        """
        # Build GFF/GBK-compatible SeqRecords for all contigs
        prodigal_records: List[SeqRecord] = []
        for seq in self.contigs:
            record = SeqRecord(seq.seq, id=seq.id, name=getattr(seq, "name", seq.id), description=getattr(seq, "description", ""))
            record.annotations["molecule_type"] = "DNA"
            record.features = []  # type: ignore[attr-defined]
            prodigal_records.append(record)

        # Train GeneFinder (prefer the user-provided assembly; fallback to reference GBFF if available)
        training_seqs: List[str] = []
        # Prefer training on the assembly contigs (best fit for sample-specific coding statistics)
        if self.contigs:
            try:
                training_seqs = [str(rec.seq) for rec in self.contigs]
            except Exception:
                training_seqs = []
        # Fallback to reference GBFF only if we couldn't extract assembly sequences
        if not training_seqs and self.reference_gbff and isinstance(self.reference_gbff, str) and os.path.isfile(self.reference_gbff):
            try:
                training_seqs = [str(rec.seq) for rec in SeqIO.parse(self.reference_gbff, "genbank")]
            except Exception:
                training_seqs = []
        gene_finder = pyrodigal.GeneFinder()
        if training_seqs:
            gene_finder.train("".join(training_seqs))

        log("Finding ORFs in assembly contigs...")
        predicted_features: List[Tuple[SeqFeature, str]] = []
        for seq_record in self.contigs:
            genes = gene_finder.find_genes(str(seq_record.seq))
            # locate corresponding prodigal record
            record = next((r for r in prodigal_records if r.id == seq_record.id), None)
            if record is None:
                log(f"Warning: No matching record found for contig {seq_record.id}")
                continue

            for gene in genes:
                # Adjust coordinates for forward strand to be 0-based half-open
                if gene.strand == 1:
                    location = FeatureLocation(gene.begin - 1, gene.end, strand=gene.strand)
                else:
                    location = FeatureLocation(gene.begin, gene.end, strand=gene.strand)

                translation = gene.translate()
                qualifiers = {
                    "translation": [translation],
                    "ID": [f"{seq_record.id}_cds_{gene.begin}_{gene.end}"],
                }
                feature = SeqFeature(location=location, type="CDS", qualifiers=qualifiers)
                record.features.append(feature)
                predicted_features.append((feature, translation))

        log(f"Found {len(predicted_features)} predicted CDS features.")
        return prodigal_records, predicted_features

    def _direct_match_protein(self, feature: SeqFeature, seq: str) -> Tuple[str, object]:
        """Return (status, annotated_feature or (feature, seq))."""
        norm = seq.replace("*", "").strip()
        if norm in self.exact_ref_lookup:
            ref_id = self.exact_ref_lookup[norm]
            ref_feature = self.protein_feature_map.get(ref_id)
            if ref_feature:
                for key in ["locus_tag", "gene", "product", "note"]:
                    if key in ref_feature.qualifiers:
                        feature.qualifiers[key] = ref_feature.qualifiers[key]
                return ("annotated", feature)
        return ("unmatched", (feature, seq))

    def _hmm_search_protein(
        self, pred_feature: SeqFeature, pred_seq: str
    ) -> Tuple[Optional[SeqFeature], List[SeqRecord], Optional[str]]:
        query = easel.TextSequence(name=b"query", sequence=pred_seq)
        digital_query = query.digitize(self.prot_alphabet)
        results = hmmer.phmmer(digital_query, list(self.digital_refs.values()))

        annotated: Optional[SeqFeature] = None
        variants: List[SeqRecord] = []
        matched: Optional[str] = None

        for hit_list in results:
            if len(hit_list) > 0:
                top_hit = hit_list[0]
                domain = top_hit.best_domain
                if domain.i_evalue > self.e_value:
                    return (None, [], None)
                ref_locus = top_hit.name.decode()
                ref_feature = self.protein_feature_map.get(ref_locus)
                if ref_feature:
                    for key in ["locus_tag", "gene", "product", "note"]:
                        if key in ref_feature.qualifiers:
                            pred_feature.qualifiers[key] = ref_feature.qualifiers[key]

                    ref_prot = ref_feature.qualifiers.get("translation", [""])[0]
                    if str(pred_seq).rstrip("*").strip() != str(ref_prot).rstrip("*").strip():
                        variants = [
                            SeqRecord(Seq(str(pred_seq).rstrip("*")), id=f"prodigal_{ref_locus}", description=""),
                            SeqRecord(Seq(str(ref_prot).rstrip("*")), id=f"reference_{ref_locus}", description=""),
                        ]
                annotated = pred_feature
                matched = ref_locus
                break
        return (annotated, variants, matched)

    def _hmm_search_rna(self, rna_id: str, digital_query: easel.DigitalSequence, contigs: List[SeqRecord]) -> List[Tuple[str, str, int, int, int]]:
        results = hmmer.nhmmer(digital_query, [easel.TextSequence(name=rec.id.encode(), sequence=str(rec.seq)).digitize(self.dna_alphabet) for rec in contigs])
        hits: List[Tuple[str, str, int, int, int]] = []
        for hit_list in results:
            if len(hit_list) > 0:
                top_hit = hit_list[0]
                contig_name = top_hit.name.decode()
                domain = top_hit.best_domain
                start = domain.env_from
                end = domain.env_to
                strand = 1 if start < end else -1
                hits.append((rna_id, contig_name, start, end, strand))
        return hits

    def annotate(self) -> List[SeqRecord]:
        """Run the full pipeline and write the annotated GBFF to disk.
        Returns the list of annotated SeqRecords.
        """
        prodigal_records, predicted_features = self._find_orfs()

        # 1) Direct matches
        annotated_features: List[SeqFeature] = []
        unmatched: Dict[str, Tuple[SeqFeature, str]] = {}
        variant_records: List[SeqRecord] = []

        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = [executor.submit(self._direct_match_protein, ftr, seq) for (ftr, seq) in predicted_features]
            with tqdm(total=len(futures), desc="[orfmatch] Checking for direct sequence matches", unit="cds") as pbar:
                for future in as_completed(futures):
                    status, result = future.result()
                    if status == "annotated":
                        annotated_features.append(result)  # type: ignore[arg-type]
                    else:
                        feature, seq = result  # type: ignore[assignment]
                        unmatched[feature.qualifiers["ID"][0]] = (feature, seq)
                    pbar.update(1)
        log(f"Found {len(annotated_features)} direct sequence matches")

        # 2) HMMER for proteins
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = [executor.submit(self._hmm_search_protein, ftr, seq) for (ftr, seq) in unmatched.values()]
            with tqdm(total=len(futures), desc="[orfmatch] Annotating unmatched CDSs using pyhmmer", unit="cds") as pbar:
                for future in as_completed(futures):
                    annotated, variants, matched = future.result()
                    if annotated:
                        annotated_features.append(annotated)
                    if variants:
                        variant_records.extend(variants)
                    if matched and annotated:
                        # Remove by the annotated feature ID
                        unmatched.pop(annotated.qualifiers["ID"][0], None)
                    pbar.update(1)

        # 3) Label remaining unmatched as hypothetical
        log("Labelling remaining unmatched CDSs as 'hypothetical protein'...")
        for feature_id, (feature, _seq) in unmatched.items():
            if feature.type == "CDS":
                if "product" not in feature.qualifiers:
                    feature.qualifiers["product"] = ["hypothetical protein"]
                feature.qualifiers["note"] = ["No match found during annotation"]
                annotated_features.append(feature)

        # 4) HMMER for RNAs
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = [
                executor.submit(self._hmm_search_rna, rna_id, digitals[0], self.contigs)
                for rna_id, digitals in self.digital_rnas.items()
                if digitals
            ]
            with tqdm(total=len(futures), desc="[orfmatch] Annotating RNAs using pyhmmer", unit="rna") as pbar:
                for future in as_completed(futures):
                    hits = future.result()
                    if hits:
                        rna_id, contig_name, start, end, strand = hits[0]
                        start0 = min(start, end) - 1
                        end0 = max(start, end)
                        location = FeatureLocation(start0, end0, strand=strand)
                        rna_feature_tpl = self.rna_feature_map[rna_id][0]
                        qualifiers = dict(rna_feature_tpl.qualifiers)
                        feature_type = rna_feature_tpl.type
                        feature = SeqFeature(location=location, type=feature_type, qualifiers=qualifiers)
                        for record in prodigal_records:
                            if record.id == contig_name:
                                record.features.append(feature)
                    pbar.update(1)

        # 5) Merge features and write output
        for record in prodigal_records:
            original_rna_features = [f for f in record.features if f.type in {"tRNA", "rRNA", "ncRNA"}]
            matched_cds_features = [
                f for f in annotated_features
                if hasattr(f, "location") and f.location is not None and record.id in f.qualifiers.get("ID", [""])[0]
            ]
            record.features = original_rna_features + matched_cds_features
            record.features.sort(key=lambda f: min(int(f.location.start), int(f.location.end)))

        with open(self.annotated_gbff, "w") as out_handle:
            SeqIO.write(prodigal_records, out_handle, "genbank")

        # 6) Variants (optional)
        if self.show_variants and variant_records:
            SeqIO.write(variant_records, self.variants_fasta, "fasta")
            log(f"  Variants found: {len(variant_records) // 2}")
            log(f"[✓] Variants saved to: {self.variants_fasta}")
            with open("variant_alignments.txt", "w") as aln_out:
                aligner = PairwiseAligner()
                aligner.mode = "global"
                for i in tqdm(range(0, len(variant_records), 2), desc="Aligning variants", unit="aln"):
                    prodigal_record = variant_records[i]
                    reference_record = variant_records[i + 1]
                    prodigal_seq = str(prodigal_record.seq).rstrip("*")
                    reference_seq = str(reference_record.seq).rstrip("*")
                    alignment = aligner.align(prodigal_seq, reference_seq)[0]
                    aln_out.write(f"Alignment of {prodigal_record.id} and {reference_record.id}: \n")
                    aln_out.write(str(alignment) + "\n")
            log(f"[✓] Writing pairwise alignments of variants to alignments.txt")

        # 7) Summary
        total_identified_rnas = sum(
            1 for record in prodigal_records for feature in record.features if feature.type in {"tRNA", "rRNA", "ncRNA"}
        )
        log("[Summary]")
        log(f"  Total reference proteins: {len(self.reference_proteins)}")
        log(f"  Total predicted proteins: {len(predicted_features)}")
        log(f"  Total reference RNAs: {len(self.reference_rnas)}")
        log(f"  Total identified RNAs: {total_identified_rnas}\n")
        log(f"[✓] Annotated GBFF written to: {self.annotated_gbff}\n")

        return prodigal_records
