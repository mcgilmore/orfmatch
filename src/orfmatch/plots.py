from pycirclize.parser import Genbank
from pycirclize import Circos
# TODO: Different genome plots using pyGenomeViz


class Circle:
    def __init__(self, reference, assembly):
        
        self.reference = reference
        self.assembly = assembly

        self.ref_gbk = Genbank(self.reference)
        self.asm_gbk = Genbank(self.assembly)

    def plot(self, dpi):
        circos = Circos(
            sectors=dict(**self.ref_gbk.get_seqid2size(), **dict(reversed(list(self.asm_gbk.get_seqid2size().items())))),
            start=-358,
            end=2,
            space=1,
            sector2clockwise={seqid: False for seqid in self.asm_gbk.get_seqid2size().keys()}, # TODO: Figure out what this actually does
        ) 

        ref_features = self.ref_gbk.get_seqid2features()
        assembly_features = self.asm_gbk.get_seqid2features()

        for sector in circos.sectors:
            cds_track = sector.add_track((95, 100))
            cds_track.axis(fc="none", ec="none")
            ref_sector_features = ref_features[sector.name]
            assembly_sector_features = assembly_features[sector.name]
            for feature in ref_sector_features:
                if feature.location.strand == 1:
                    cds_track.genomic_features(feature, plotstyle="arrow", r_lim=(95, 100), fc="salmon")
                else:
                    cds_track.genomic_features(feature, plotstyle="arrow", r_lim=(95, 100), fc="skyblue")
            for feature in assembly_sector_features:
                if feature.location.strand == 1:
                    cds_track.genomic_features(feature, plotstyle="arrow", r_lim=(95, 100), fc="salmon")
                else:
                    cds_track.genomic_features(feature, plotstyle="arrow", r_lim=(95, 100), fc="skyblue")
            
        
        # Build a lookup: locus_tag -> (sector_name, feature)
        ref_locus_map = {}
        for seqid, features in ref_features.items():
            for f in features:
                locus = f.qualifiers.get("locus_tag", [""])[0]
                ref_locus_map[locus] = (seqid, f)

        asm_locus_map = {}
        for seqid, features in assembly_features.items():
            for f in features:
                locus = f.qualifiers.get("locus_tag", [""])[0]
                asm_locus_map[locus] = (seqid, f)

        # Now find and link matches
        for locus in set(ref_locus_map.keys()) & set(asm_locus_map.keys()):
            ref_seqid, ref_feat = ref_locus_map[locus]
            asm_seqid, asm_feat = asm_locus_map[locus]
            circos.link(
                (ref_seqid, int(ref_feat.location.start), int(ref_feat.location.end)),
                (asm_seqid, int(asm_feat.location.start), int(asm_feat.location.end)),
                color="blue"
            )

        circos.savefig("plot.png", dpi=dpi)