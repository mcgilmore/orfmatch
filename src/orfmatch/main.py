import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pyrodigal
from collections import defaultdict
from orfmatch.utils import log
from orfmatch.annotation import Annotator
import os
from importlib.metadata import version


def main():
    log(f"version {version('orfmatch')}")

    parser = argparse.ArgumentParser(
        description="Transfer feature annotations from a reference genome to a de novo assembled one.")
    parser.add_argument("-v", "--variants", action="store_true",
                        help="Output protein sequences which differ from reference genome")
    parser.add_argument("-e", "--e-value", type=float, default=1e-25,
                        help="E-value threshold for accepting phmmer matches (default: 1e-25)")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of threads for parallel steps (default: 8)")
    parser.add_argument("-c", "--circle", action="store_true",
                        help="Output circular plot of reference against assembly annotations to SVG")
    parser.add_argument("-l", "--line", action="store_true",
                        help="Output linear plot of reference against assembly annotations to SVG")
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA assembly")
    parser.add_argument("-r", "--reference", required=True,
                        help="Reference GBFF file with annotations")
    parser.add_argument("-o", "--output", required=True, 
                        help="Output GFF file with transferred annotations")
    parser.add_argument("-u", "--unique", action="store_true", # TODO: Output unique features
                        help="List features which are present only in the reference or query")

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        parser.error(f"Input FASTA file '{args.input}' not found.")
    if not os.path.isfile(args.reference):
        parser.error(f"Reference GBFF file '{args.reference}' not found.")

    assembly_fasta = args.input
    reference_gbff = args.reference
    annotated_gbff = args.output
    show_variants = args.variants

    # Sanitise input # TODO: Protein fasta mode
    output_base, ext = os.path.splitext(annotated_gbff)
    if ext != ".gbff":
        annotated_gbff = output_base + ".gbff"

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

    log("Parsing reference sequence...")
    rna_id_counter = defaultdict(int)
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

                rna_id_counter[locus] += 1
                if rna_id_counter[locus] > 1:
                    locus = f"{locus}_{rna_id_counter[locus]}"

                reference_rnas.append(
                    SeqRecord(rna_seq, id=locus, description=feature.type))
                if locus not in rna_feature_map:
                    rna_feature_map[locus] = []
                rna_feature_map[locus].append(feature)
    exact_ref_lookup.update({str(p.seq): p.id for p in reference_proteins})
    log(f"Found {len(reference_proteins)} ORFs and {len(reference_rnas)} RNAs.")

    annotator = Annotator()
    annotater.annotate()

    if args.detect_novel:
        from orfmatch.alignment import Aligner
        
        log("Finding regions not present in reference...")
        novel = Aligner().find_novel_regions(contigs, reference_gbff, min_length=100)

        for record in prodigal_records:
            for start, end in novel.get(record.id, []):
                location = FeatureLocation(start, end, strand=1)
                feature = SeqFeature(
                    location=location,
                    type="novel_region",
                    qualifiers={"note": ["No alignment to reference genome"]}
                )
                record.features.append(feature)
        log("[✓] Novel regions added to feature list")

    if args.circle:
        log(f"[✓] Plotting circular comparison and saving to {output_base}_circle_plot.svg")
        circle = Circle(reference_gbff, annotated_gbff)
        circle.plot(f"{output_base}_circle_plot")
    if args.line:
        log(f"[✓] Plotting linear comparison and saving to {output_base}_line_plot.svg")
        line = Line(reference_gbff, annotated_gbff)
        line.plot(f"{output_base}_line_plot")


if __name__ == "__main__":
    main()
