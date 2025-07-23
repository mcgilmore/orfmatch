class Aligner:
    import mappy as mp
    from Bio import SeqRecord, Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    assembly_fasta = ""
    reference_fasta = ""

    def find_novel_regions(query_records, reference_fasta, min_length=100):
        aligner = mp.Aligner(reference_fasta, preset="asm5")
        if not aligner:
            raise RuntimeError("Failed to load reference for alignment.")

        novel_regions = {}  # dict[contig_id] = list[(start, end)]

        for record in query_records:
            seq = str(record.seq)
            covered = []

            for hit in aligner.map(seq):
                covered.append((hit.q_st, hit.q_en))

            covered.sort()
            gaps = []
            pos = 0
            for start, end in covered:
                if start > pos and (start - pos) >= min_length:
                    gaps.append((pos, start))
                pos = max(pos, end)
            if pos < len(seq) and (len(seq) - pos) >= min_length:
                gaps.append((pos, len(seq)))

            novel_regions[record.id] = gaps

        return novel_regions
