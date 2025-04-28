from pycirclize.parser import Genbank
from pycirclize import Circos
import seaborn as sns
# TODO: Different genome plots using pyGenomeViz


class Circle:
    def __init__(self, reference, assembly):
        self.ref_gbk = Genbank(reference)
        self.asm_gbk = Genbank(assembly)

    def plot(self):
        # Rename contigs to avoid problems with identically named ones
        sectors = {}
        sectors.update(self.ref_gbk.get_seqid2size())
        sectors.update({f"ASM_{k}": v for k, v in reversed(list(self.asm_gbk.get_seqid2size().items()))})
        
        circos = Circos(
            sectors=sectors,
            start=-358,
            end=2,
            space=3,
            sector2clockwise={
                seqid: False for seqid in self.asm_gbk.get_seqid2size().keys()},
        )

        ref_features = self.ref_gbk.get_seqid2features()
        assembly_features = self.asm_gbk.get_seqid2features()

        for sector in circos.sectors:
            label = sector.name
            if label.startswith("ASM_"):
                label = label.replace("ASM_", "")
            sector.text(
                label,
                r=61,
                size=6
            )
            cds_track = sector.add_track((59.8, 60.2))
            cds_track.axis(fc="black", ec="none")
            
            if sector.name.startswith("ASM_"):
                asm_name = sector.name.replace("ASM_", "")
                ref_sector_features = []
                assembly_sector_features = assembly_features.get(asm_name, [])
            else:
                ref_sector_features = ref_features.get(sector.name, [])
                assembly_sector_features = []
            
            for feature in ref_sector_features:
                if feature.location.strand == 1:
                    cds_track.genomic_features(
                        feature, plotstyle="arrow", r_lim=(59.8, 60.2), fc="salmon")
                else:
                    cds_track.genomic_features(
                        feature, plotstyle="arrow", r_lim=(59.8, 60.2), fc="skyblue")
            for feature in assembly_sector_features:
                if feature.location.strand == 1:
                    cds_track.genomic_features(
                        feature, plotstyle="arrow", r_lim=(59.8, 60.2), fc="salmon")
                else:
                    cds_track.genomic_features(
                        feature, plotstyle="arrow", r_lim=(59.8, 60.2), fc="skyblue")

        # Build a lookup: locus_tag -> (sector_name, feature)
        ref_locus_map = {}
        for seqid, features in ref_features.items():
            for f in features:
                locus = f.qualifiers.get("locus_tag", [""])[0]
                ref_locus_map[locus] = (seqid, f)

        asm_locus_map = {}
        for seqid, features in assembly_features.items():
            asm_seqid = f"ASM_{seqid}"  # <- Add prefix here
            for f in features:
                locus = f.qualifiers.get("locus_tag", [""])[0]
                asm_locus_map[locus] = (asm_seqid, f)

        # Assign a color per reference contig using seaborn
        ref_contig_colors = {}
        ref_contigs = sorted(set(seqid for seqid, _ in ref_locus_map.values()))
        palette = sns.color_palette("muted", n_colors=len(ref_contigs))

        for idx, seqid in enumerate(ref_contigs):
            ref_contig_colors[seqid] = palette[idx]

        # Now create links using the color per reference contig
        for locus in set(ref_locus_map.keys()) & set(asm_locus_map.keys()):
            ref_seqid, ref_feat = ref_locus_map[locus]
            asm_seqid, asm_feat = asm_locus_map[locus]
            color = ref_contig_colors.get(ref_seqid, "lightgray")
            circos.link(
                (ref_seqid, int(ref_feat.location.start), int(ref_feat.location.end)),
                (asm_seqid, int(asm_feat.location.start), int(asm_feat.location.end)),
                color=color
            )

        # Add Reference and Assembly side labels
        circos.text(
            "Reference",
            r=72,
            deg=90,
            size=10,
            orientation="vertical",
            adjust_rotation=True
        )
        circos.text(
            "Assembly",
            r=72,
            deg=-90,
            size=10,
            orientation="vertical",
            adjust_rotation=True
        )

        circos.savefig("circle.svg")
