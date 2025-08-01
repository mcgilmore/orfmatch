from __future__ import annotations

from typing import Dict, List, Tuple
from collections import defaultdict
import os
import tempfile

import mappy as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


class Aligner:
    """
    Align query sequences against a reference genome (from GBFF) using minimap2
    and emit **DNA FASTA files** of differences (no GBFF annotation).

    Outputs (created only if non-empty):
      - <prefix>_query_novel.fa         : query segments not covered by any alignment (len ≥ min_length_novel)
      - <prefix>_query_insertions.fa    : query insertion blocks (CIGAR I; len ≥ min_length_variant)
      - <prefix>_reference_deletions.fa : reference deletion/skip blocks (CIGAR D/N; len ≥ min_length_variant)
      - <prefix>_query_mismatches.fa    : query mismatch blocks (CIGAR X; len ≥ min_length_variant)
      - <prefix>_reference_mismatches.fa: reference mismatch blocks (CIGAR X; len ≥ min_length_variant)
    """

    def __init__(self, preset: str = "asm5") -> None:
        self.preset = preset

    # ------------------------ helpers ------------------------
    @staticmethod
    def _write_temp_fasta_from_gbff(reference_gbff: str) -> str:
        """Extract sequences from a GBFF to a temporary FASTA file and return the path."""
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".fa")
        tmp_path = tmp.name
        tmp.close()
        with open(tmp_path, "w") as out_fa:
            for rec in SeqIO.parse(reference_gbff, "genbank"):
                out_fa.write(f">{rec.id}\n{str(rec.seq)}\n")
        return tmp_path

    @staticmethod
    def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        if not intervals:
            return []
        intervals.sort()
        merged = [intervals[0]]
        for s, e in intervals[1:]:
            ps, pe = merged[-1]
            if s <= pe:
                merged[-1] = (ps, max(pe, e))
            else:
                merged.append((s, e))
        return merged

    @staticmethod
    def _compute_gaps(length: int, covered: List[Tuple[int, int]], min_len: int) -> List[Tuple[int, int]]:
        """Return uncovered (start, end) intervals of the query given covered spans."""
        gaps: List[Tuple[int, int]] = []
        if length <= 0:
            return gaps
        merged = Aligner._merge_intervals(covered)
        pos = 0
        for s, e in merged:
            if s > pos and (s - pos) >= min_len:
                gaps.append((pos, s))
            pos = max(pos, e)
        if pos < length and (length - pos) >= min_len:
            gaps.append((pos, length))
        return gaps

    @staticmethod
    def _revcomp(seq: str) -> str:
        tbl = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return seq.translate(tbl)[::-1]

    def write_difference_fastas(
        self,
        query_records: List[SeqRecord],
        reference_gbff: str,
        *,
        min_length_novel: int = 100,
        min_length_variant: int = 1,
        out_prefix: str = "differences",
        rc_reference_minus_strand: bool = False,
    ) -> Dict[str, str]:
        """
        Align queries to reference and write DNA FASTA files for novel regions,
        insertions, deletions/skips, and mismatches. Returns a dict of kind->path
        for files that were written.
        """
        # Load reference sequences from GBFF
        ref_records: List[SeqRecord] = list(SeqIO.parse(reference_gbff, "genbank"))
        ref_seq_by_id: Dict[str, str] = {r.id: str(r.seq) for r in ref_records}

        # Prepare containers
        query_novel: List[Tuple[str, str]] = []
        query_insertions: List[Tuple[str, str]] = []
        query_mismatches: List[Tuple[str, str]] = []
        ref_deletions: List[Tuple[str, str]] = []
        ref_mismatches: List[Tuple[str, str]] = []

        # Helper to write FASTA
        def write_fasta(path: str, entries: List[Tuple[str, str]]) -> None:
            with open(path, "w") as fh:
                for header, seq in entries:
                    fh.write(">" + header + "\n")
                    fh.write(seq + "\n")

        # Align
        ref_fa = self._write_temp_fasta_from_gbff(reference_gbff)
        try:
            aligner = mp.Aligner(ref_fa, preset=self.preset)
            if not aligner:
                raise RuntimeError("Failed to load reference for alignment.")

            for qrec in query_records:
                qseq = str(qrec.seq)
                covered: List[Tuple[int, int]] = []

                for hit in aligner.map(qseq):
                    covered.append((hit.q_st, hit.q_en))

                    qpos = hit.q_st
                    rpos = hit.r_st
                    ref_id = hit.ctg
                    strand = hit.strand  # 1 or -1
                    rseq = ref_seq_by_id.get(ref_id)
                    if rseq is None:
                        continue

                    cigar = getattr(hit, "cigar", None)
                    if cigar is None:
                        continue

                    for op, length in cigar:
                        if length <= 0:
                            continue
                        if op in (0, 7):  # M or =
                            qpos += length
                            rpos = rpos + length if strand >= 0 else rpos - length
                        elif op == 8:  # X mismatch block
                            if length >= min_length_variant:
                                # query side
                                q_start, q_end = qpos, qpos + length
                                q_seq = qseq[q_start:q_end]
                                q_hdr = f"{qrec.id}|mismatch|start={q_start}|end={q_end}|ref={ref_id}|strand={strand}"
                                query_mismatches.append((q_hdr, q_seq))
                                # reference side
                                if strand >= 0:
                                    r_start, r_end = rpos, rpos + length
                                else:
                                    r_start, r_end = rpos - length, rpos
                                r_seq = rseq[r_start:r_end]
                                if strand < 0 and rc_reference_minus_strand:
                                    r_seq = self._revcomp(r_seq)
                                r_hdr = f"{ref_id}|mismatch|start={r_start}|end={r_end}|query={qrec.id}|strand={strand}"
                                ref_mismatches.append((r_hdr, r_seq))
                            qpos += length
                            rpos = rpos + length if strand >= 0 else rpos - length
                        elif op == 1:  # I insertion w.r.t reference (extra bases in query)
                            if length >= min_length_variant:
                                q_start, q_end = qpos, qpos + length
                                q_seq = qseq[q_start:q_end]
                                q_hdr = f"{qrec.id}|insertion|start={q_start}|end={q_end}|ref={ref_id}|strand={strand}"
                                query_insertions.append((q_hdr, q_seq))
                            qpos += length
                        elif op == 2 or op == 3:  # D or N on reference
                            if length >= min_length_variant:
                                if strand >= 0:
                                    r_start, r_end = rpos, rpos + length
                                else:
                                    r_start, r_end = rpos - length, rpos
                                r_seq = rseq[r_start:r_end]
                                if strand < 0 and rc_reference_minus_strand:
                                    r_seq = self._revcomp(r_seq)
                                r_hdr = f"{ref_id}|deletion|start={r_start}|end={r_end}|query={qrec.id}|strand={strand}"
                                ref_deletions.append((r_hdr, r_seq))
                            rpos = rpos + length if strand >= 0 else rpos - length
                        elif op == 4:  # soft-clip consumes query only
                            qpos += length
                        elif op == 5:  # hard-clip consumes neither sequence here
                            continue
                        else:
                            # fallback: advance if it consumes query
                            qpos += length

                # Novel regions on the query (uncovered by any alignment)
                for start, end in self._compute_gaps(len(qseq), covered, min_length_novel):
                    q_hdr = f"{qrec.id}|novel|start={start}|end={end}"
                    query_novel.append((q_hdr, qseq[start:end]))

            # Write FASTAs
            written: Dict[str, str] = {}
            if query_novel:
                p = f"{out_prefix}_query_novel.fa"
                write_fasta(p, query_novel)
                written["query_novel"] = p
            if query_insertions:
                p = f"{out_prefix}_query_insertions.fa"
                write_fasta(p, query_insertions)
                written["query_insertions"] = p
            if ref_deletions:
                p = f"{out_prefix}_reference_deletions.fa"
                write_fasta(p, ref_deletions)
                written["reference_deletions"] = p
            if query_mismatches:
                p = f"{out_prefix}_query_mismatches.fa"
                write_fasta(p, query_mismatches)
                written["query_mismatches"] = p
            if ref_mismatches:
                p = f"{out_prefix}_reference_mismatches.fa"
                write_fasta(p, ref_mismatches)
                written["reference_mismatches"] = p

            return written
        finally:
            try:
                os.unlink(ref_fa)
            except Exception:
                pass
