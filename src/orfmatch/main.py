import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Align import PairwiseAligner
import pyrodigal
from pyhmmer import easel, hmmer
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

# TODO: Filtering of hits by E value to get rid of those which are unlikely to be true

def direct_match(predicted):
    feature, seq = predicted
    for ref_seq, ref_id in exact_ref_lookup.items():
        if seq.replace("*", "").strip() == ref_seq.replace("*", "").strip():
            ref_feature = protein_feature_map.get(ref_id)
            if ref_feature:
                for key in ["locus_tag", "gene", "product", "note"]:
                    if key in ref_feature.qualifiers:
                        feature.qualifiers[key] = ref_feature.qualifiers[key]
                return ("annotated", feature)
    return ("unmatched", (feature, seq))


def search_and_annotate(pred_feature, pred_seq, refs, feature_map, alphabet):
    query = easel.TextSequence(name=b"query", sequence=pred_seq)
    digital_query = query.digitize(alphabet)
    results = hmmer.phmmer(digital_query, list(refs.values()))

    annotated = None
    variants = []
    matched = None
    for hit_list in results:
        if len(hit_list) > 0:
            ref_locus = hit_list[0].name.decode()
            ref_feature = feature_map.get(ref_locus)
            if ref_feature:
                for key in ["locus_tag", "gene", "product", "note"]:
                    if key in ref_feature.qualifiers:
                        pred_feature.qualifiers[key] = ref_feature.qualifiers[key]

                ref_prot = ref_feature.qualifiers["translation"][0]
                # Check for exact match, and if not, add to variants
                if str(pred_seq).rstrip("*").strip() != ref_prot.rstrip("*").strip():
                    variants = [
                        SeqRecord(Seq(pred_seq.rstrip("*")),
                                  id=f"prodigal_{ref_locus}", description=""),
                        SeqRecord(Seq(ref_prot.rstrip("*")),
                                  id=f"reference_{ref_locus}", description="")
                    ]
            annotated = pred_feature
            matched = ref_locus
            break
    return (annotated, variants, matched)


def search_rna_annotation(rna_id_query, contigs, dna_alphabet):
    rna_id, query = rna_id_query
    results = hmmer.nhmmer(query, [easel.TextSequence(name=rec.id.encode(
    ), sequence=str(rec.seq)).digitize(dna_alphabet) for rec in contigs])
    hits = []
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


exact_ref_lookup = {}


def main():
    parser = argparse.ArgumentParser(
        description="Annotate Prodigal-predicted CDSs using a reference GBFF and pyhmmer.")
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA assembly")
    parser.add_argument("-r", "--reference", required=True,
                        help="Reference GBFF file with annotations")
    parser.add_argument("-o", "--output", required=True,
                        help="Output GFF file with transferred annotations")
    parser.add_argument("-v", "--variants", action="store_true",
                        help="Output protein sequences which differ from reference genome")
    args = parser.parse_args()

    assembly_fasta = args.input
    reference_gbff = args.reference
    annotated_gbff = args.output
    show_variants = args.variants

    # Load all contigs from the assembly FASTA
    contigs = list(SeqIO.parse(assembly_fasta, "fasta"))

    # Outputs
    variants_fasta = "variants.fasta"

    # Step 1: Extract reference protein sequences
    reference_proteins = []
    global protein_feature_map
    protein_feature_map = {}

    reference_rnas = []
    rna_feature_map = {}

    print("\nParsing reference sequence...", end=" ")
    for record in SeqIO.parse(reference_gbff, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                prot_seq = feature.qualifiers["translation"][0]
                locus = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                protein = SeqRecord(Seq(prot_seq), id=locus, description="")
                reference_proteins.append(protein)
                protein_feature_map[locus] = feature
            elif feature.type in {"tRNA", "rRNA", "ncRNA"}:
                rna_seq = feature.extract(record.seq)
                locus = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                reference_rnas.append(
                    SeqRecord(rna_seq, id=locus, description=feature.type))
                if locus not in rna_feature_map:
                    rna_feature_map[locus] = []
                rna_feature_map[locus].append(feature)
    exact_ref_lookup.update({str(p.seq): p.id for p in reference_proteins})
    print(
        f"Found {len(reference_proteins)} ORFs and {len(reference_rnas)} RNAs.")

    # Step 2: Convert reference proteins to digital sequences for phmmer
    alphabet = easel.Alphabet.amino()
    digital_refs = {
        rec.id: easel.TextSequence(
            name=rec.id.encode(), sequence=str(rec.seq)).digitize(alphabet)
        for rec in reference_proteins
    }

    # Convert reference RNAs to digital sequences for nhmmer
    dna_alphabet = easel.Alphabet.dna()
    digital_rnas = defaultdict(list)
    for rec in reference_rnas:
        digital = easel.TextSequence(name=rec.id.encode(), sequence=str(rec.seq)).digitize(dna_alphabet)
        digital_rnas[rec.id].append(digital)

    # Generate GFF-compatible records with feature lists
    prodigal_records = []
    for seq in contigs:
        record = SeqRecord(seq.seq, id=seq.id, name=seq.name,
                           description=seq.description)
        record.annotations["molecule_type"] = "DNA"
        prodigal_records.append(record)

    for record in prodigal_records:
        record.features = []

    # Extract nucleotide sequences from reference GBFF for training
    training_seqs = []
    for record in SeqIO.parse(reference_gbff, "genbank"):
        training_seqs.append(str(record.seq))

    # Train on concatenated genome if multiple records
    gene_finder = pyrodigal.GeneFinder()
    gene_finder.train("".join(training_seqs))

    predicted_features = []
    variant_records = []

    print("Finding ORFs in assembly contigs...", end=" ")
    for seq_record in contigs:
        genes = gene_finder.find_genes(str(seq_record.seq))
        record = next(
            (r for r in prodigal_records if r.id == seq_record.id), None)
        if not record:
            print(
                f"Warning: No matching record found for contig {seq_record.id}")
            continue

        for gene in genes:
            if gene.strand == 1:
                # Forward strand CDSs need to be adjusted
                location = FeatureLocation(
                    gene.begin - 1, gene.end, strand=gene.strand)
            else:
                location = FeatureLocation(
                    gene.begin, gene.end, strand=gene.strand)

            qualifiers = {
                "translation": [gene.translate()],
                # TODO: Probably remove this ID
                "ID": [f"{seq_record.id}_cds_{gene.begin}_{gene.end}"]
            }
            feature = SeqFeature(
                location=location, type="CDS", qualifiers=qualifiers)
            record.features.append(feature)
            predicted_features.append((feature, gene.translate()))
    print(f"Found {len(predicted_features)}.")

    # Step 4: Directly match genes with their annotation for identical sequences
    annotated_features = []
    unmatched = {}
    variant_records = []

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(direct_match, p)
                   for p in predicted_features]

        with tqdm(total=len(futures), desc="Checking for direct sequence matches", unit="cds") as pbar:
            for future in as_completed(futures):
                status, result = future.result()
                if status == "annotated":
                    annotated_features.append(result)
                else:
                    feature, seq = result
                    unmatched[feature.qualifiers["ID"][0]] = (feature, seq)
                pbar.update(1)
    print(f"Found {len(annotated_features)} direct sequence matches")

    # Step 4.5: Search with phmmer (parallelized)
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(search_and_annotate, ftr, s, digital_refs, protein_feature_map, alphabet)
                   for ftr, (ftr, s) in unmatched.items()]

        with tqdm(total=len(futures), desc="Annotating unmatched CDSs using pyhmmer", unit="cds") as pbar:
            for future in as_completed(futures):
                annotated, variants, matched = future.result()
                if annotated:
                    annotated_features.append(annotated)
                if variants:
                    variant_records.extend(variants)
                if matched:
                    digital_refs.pop(matched, None)
                    # Also remove the matched feature from unmatched list
                    unmatched.pop(annotated.qualifiers["ID"][0], None)
                pbar.update(1)

    # Step 4.55: Label unmatched proteins as hypothetical and add them
    print("Labelling remaining unmatched CDSs as 'hypothetical protein'...")
    for feature_id, (feature, seq) in unmatched.items():
        if feature.type == "CDS":
            if "product" not in feature.qualifiers:
                feature.qualifiers["product"] = ["hypothetical protein"]
            feature.qualifiers["note"] = ["No match found during annotation"]
            annotated_features.append(feature)

    # Step 4.6: Search with nhmmer (for RNAs)
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(
            search_rna_annotation, (rna_id, digitals[0]), contigs, dna_alphabet) 
            for rna_id, digitals in digital_rnas.items()]

        with tqdm(total=len(futures), desc="Annotating RNAs using pyhmmer", unit="rna") as pbar:
            for future in as_completed(futures):
                hits = future.result()
                if hits:
                    # Take only the best hit
                    rna_id, contig_name, start, end, strand = hits[0]
                    start0 = min(start, end) - 1
                    end0 = max(start, end)
                    location = FeatureLocation(start0, end0, strand=strand)
                    rna_feature = rna_feature_map[rna_id][0]
                    feature_type = rna_feature.type
                    qualifiers = rna_feature.qualifiers.copy()
                    feature = SeqFeature(
                        location=location, type=feature_type, qualifiers=qualifiers)

                    for record in prodigal_records:
                        if record.id == contig_name:
                            record.features.append(feature)
                pbar.update(1)

    # Step 5: Output annotated GBFF
    # Output all contigs with their annotated features
    for record in prodigal_records:
        original_rna_features = [f for f in record.features if f.type in {"tRNA", "rRNA", "ncRNA"}]
        # Add new annotated CDS features
        matched_cds_features = [f for f in annotated_features if hasattr(f, "location") and hasattr(
            record, "id") and f.location is not None and record.id in f.qualifiers.get("ID", [""])[0]]

        # Merge and sort
        record.features = original_rna_features + matched_cds_features
        record.features.sort(key=lambda f: min(int(f.location.start), int(f.location.end)))
    with open(annotated_gbff, "w") as out_handle:
        SeqIO.write(prodigal_records, out_handle, "genbank")

    # Step 6: Print summary
    print("\n[Summary]")
    print(f"  Total reference proteins: {len(reference_proteins)}")
    print(f"  Total predicted proteins: {len(predicted_features)}")
    print(
        f"  Matched annotations: {len(annotated_features) - len(unmatched)}\n")
    print(f"  Total reference RNAs: {len(reference_rnas)}")
    total_identified_rnas = sum(
        1 for record in prodigal_records for feature in record.features if feature.type in {"tRNA", "rRNA", "ncRNA"}
    )
    print(f"  Total identified RNAs: {total_identified_rnas}\n")
    print(f"[✓] Annotated GBFF written to: {annotated_gbff}\n")
    if show_variants and variant_records:
        SeqIO.write(variant_records, variants_fasta, "fasta")
        print(f"  Variants found: {len(variant_records) // 2}")
        print(f"[✓] Variants saved to: {variants_fasta}")
        with open("variant_alignments.txt", "w") as aln_out:
            aligner = PairwiseAligner()
            aligner.mode = "global"
            for i in tqdm(range(0, len(variant_records), 2), desc="Aligning variants", unit="aln"):
                prodigal_record = variant_records[i]
                reference_record = variant_records[i+1]

                prodigal_seq = str(prodigal_record.seq).rstrip("*")
                reference_seq = str(reference_record.seq).rstrip("*")

                alignment = aligner.align(prodigal_seq, reference_seq)[0]
                aln_out.write(
                    f"Alignment of {prodigal_record.id} and {reference_record.id}: \n")
                aln_out.write(str(alignment) + "\n")
        print(f"[✓] Writing pairwise alignments of variants to alignments.txt")


if __name__ == "__main__":
    main()
