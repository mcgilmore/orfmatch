import argparse
from Bio import SeqIO
from orfmatch.utils import *
from orfmatch.annotation import Annotator
from orfmatch.plots import Circle, Line
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
                        help="Reference file with annotations (*.gbff or *.fasta/*.faa)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output GFF file with transferred annotations")
    parser.add_argument("-u", "--unique", action="store_true",  # TODO: Output unique features
                        help="List features which are present only in the reference or query")

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        parser.error(f"Input FASTA file '{args.input}' not found.")
    if not os.path.isfile(args.reference):
        parser.error(f"Reference file '{args.reference}' not found.")

    assembly_fasta = args.input
    reference_file = args.reference
    annotated_gbff = args.output
    show_variants = args.variants

    # Sanitise output extension
    output_base, output_ext = os.path.splitext(annotated_gbff)
    if output_ext != ".gbff":
        annotated_gbff = output_base + ".gbff"

    reference_base, reference_ext = os.path.splitext(reference_file)
    reference_ext = reference_ext.lower()

    # Load all contigs from the assembly FASTA
    contigs = list(SeqIO.parse(assembly_fasta, "fasta"))

    # Step 1: Extract reference sequences
    log(f"Parsing reference sequence: {reference_file}")
    if reference_ext == ".gbff" or reference_ext in (".gbk", ".genbank"):
        reference_proteins, reference_rnas, protein_feature_map, rna_feature_map = parse_gbff(
            reference_file)
        reference_gbff = reference_file
        log(f"Found {len(reference_proteins)} ORFs and {len(reference_rnas)} RNAs.")
    elif reference_ext in (".faa", ".fasta", ".fa"):
        reference_proteins, reference_rnas, protein_feature_map, rna_feature_map = parse_protein_fasta(
            reference_file)
        reference_gbff = ""
        log(f"Found {len(reference_proteins)} proteins.")
    else:
        parser.error(
            "Unknown reference file format. Provide *.gbff or protein FASTA (*.faa, *.fasta, *.fa).")

    # Derive a variants FASTA path next to the output
    variants_fasta = f"{output_base}_variants.faa"

    annotator = Annotator(
        contigs=contigs,
        reference_proteins=reference_proteins,
        reference_rnas=reference_rnas,
        reference_gbff=reference_gbff,
        protein_feature_map=protein_feature_map,
        rna_feature_map=rna_feature_map,
        annotated_gbff=annotated_gbff,
        variants_fasta=variants_fasta,
        threads=args.threads,
        e_value=args.e_value,
        show_variants=show_variants,
    )

    prodigal_records = annotator.annotate()

    # if args.detect_novel:
    #    from orfmatch.alignment import Aligner
    #
    #    log("Finding regions not present in reference...")
    #    novel = Aligner().find_novel_regions(contigs, reference_gbff, min_length=100)
#
    #    for record in prodigal_records:
    #        for start, end in novel.get(record.id, []):
    #            location = FeatureLocation(start, end, strand=1)
    #            feature = SeqFeature(
    #                location=location,
    #                type="novel_region",
    #                qualifiers={"note": ["No alignment to reference genome"]}
    #            )
    #            record.features.append(feature)
    #    log("[✓] Novel regions added to feature list")

    if args.circle or args.line:
        if reference_ext in (".faa", ".fasta", ".fa"):
            log("Cannot produce genome comparison plots when protein fasta given as input")
        else:
            if args.circle:
                log(
                    f"[✓] Plotting circular comparison and saving to {output_base}_circle_plot.svg")
                circle = Circle(reference_gbff, annotated_gbff)
                circle.plot(f"{output_base}_circle_plot")
            if args.line:
                log(
                    f"[✓] Plotting linear comparison and saving to {output_base}_line_plot.svg")
                line = Line(reference_gbff, annotated_gbff)
                line.plot(f"{output_base}_line_plot")


if __name__ == "__main__":
    main()
