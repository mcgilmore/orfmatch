import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF
from pyhmmer import easel, plan7, hmmer
from collections import defaultdict
import pyrodigal

def main():
    parser = argparse.ArgumentParser(description="Annotate Prodigal-predicted CDSs using a reference GBFF and pyhmmer.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA assembly")
    parser.add_argument("-r", "--reference", required=True, help="Reference GBFF file with annotations")
    parser.add_argument("-o", "--output", required=True, help="Output GFF file with transferred annotations")
    parser.add_argument("-v", "--variants", required=False, default=False, help="Output protein sequences which differ from reference genome")
    args = parser.parse_args()

    assembly_fasta = args.input
    reference_gbff = args.reference
    annotated_gff = args.output
    show_variants = args.variants

    # Load all contigs from the assembly FASTA
    contigs = list(SeqIO.parse(assembly_fasta, "fasta"))

    # Outputs
    variants_fasta = "variants.fasta"

    # Step 1: Extract reference protein sequences
    reference_proteins = []
    protein_feature_map = {}

    for record in SeqIO.parse(reference_gbff, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                prot_seq = feature.qualifiers["translation"][0]
                locus = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                protein = SeqRecord(Seq(prot_seq), id=locus, description="")
                reference_proteins.append(protein)
                protein_feature_map[locus] = feature

    # Step 2: Convert reference proteins to digital sequences for phmmer
    alphabet = easel.Alphabet.amino()
    digital_refs = []
    for rec in reference_proteins:
        text_seq = easel.TextSequence(name=rec.id.encode(), sequence=str(rec.seq))
        digital_refs.append(text_seq.digitize(alphabet))

    # Generate GFF-compatible records with feature lists
    prodigal_records = []
    for seq in contigs:
        prodigal_records.append(SeqRecord(seq.seq, id=seq.id, name=seq.name, description=seq.description))

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

    for seq_record in contigs:
        print(f"Processing contig: {seq_record.id}")
        genes = gene_finder.find_genes(str(seq_record.seq))
        record = next((r for r in prodigal_records if r.id == seq_record.id), None)
        if not record:
            print(f"Warning: No matching record found for contig {seq_record.id}")
            continue

        for gene in genes:
            print(f"Found gene: start={gene.begin}, end={gene.end}, strand={gene.strand}")
            start = gene.begin
            end = gene.end
            strand = 1 if gene.strand == 1 else -1
            location = FeatureLocation(start, end, strand=strand)
            qualifiers = {
                "translation": [gene.translate()],
                "ID": [f"{seq_record.id}_cds_{start}_{end}"]
            }
            feature = SeqFeature(location=location, type="CDS", qualifiers=qualifiers)
            record.features.append(feature)
            predicted_features.append((feature, gene.translate()))

    # Step 4: Search with phmmer
    annotated_features = []

    for pred_feature, pred_seq in predicted_features:
        query = easel.TextSequence(name=b"query", sequence=pred_seq)
        digital_query = query.digitize(alphabet)
        results = hmmer.phmmer(digital_query, digital_refs)

        top_hits = list(results) # It's a list of tophits, which are lists of hits for each gene. For each one, we take the top hit

        for hit_list in top_hits:
            if len(hit_list) > 0:
                ref_locus = hit_list[0].name
                print(f"Identified: {ref_locus}")
                ref_feature = protein_feature_map.get(ref_locus)

                if ref_feature:
                    for key in ["locus_tag", "gene", "product", "note"]:
                        if key in ref_feature.qualifiers:
                            pred_feature.qualifiers[key] = ref_feature.qualifiers[key]

                    ref_prot = ref_feature.qualifiers["translation"][0]
                    if str(pred_seq) != ref_prot:
                        variant_records.append(SeqRecord(Seq(pred_seq), id=f"prodigal_{ref_locus}", description=""))
                        variant_records.append(SeqRecord(Seq(ref_prot), id=f"reference_{ref_locus}", description=""))
            else:
                print(f"No hits found for sequence: {hit_list.query}")

        annotated_features.append(pred_feature)

    # Step 5: Output annotated GFF
    # Output all contigs with their annotated features
    for record in prodigal_records:
        record.features = [f for f in annotated_features if hasattr(f, "location") and hasattr(record, "id") and f.location is not None and record.id in f.qualifiers.get("ID", [""])[0]]
    with open(annotated_gff, "w") as out_handle:
        GFF.write(prodigal_records, out_handle)

    # Step 6: Output variants
    if show_variants and variant_records:
        SeqIO.write(variant_records, variants_fasta, "fasta")
        print(f"[✓] Variants saved to: {variants_fasta}")

    print(f"[✓] Annotated GFF written to: {annotated_gff}")


if __name__ == "__main__":
    main()