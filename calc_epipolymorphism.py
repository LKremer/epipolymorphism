#!/usr/bin/env python3


from __future__ import print_function, division
from collections import Counter, deque
import pysam
import argparse
import sys


# How to interpret CpG dinucleotides in bisulfite-converted reads:
is_meth = {"CG": 1, "TG": 0, "CA": 0}


def spans_region(read, region_start, region_end):
    """ checks whether a read fully spans a region or not """
    return read.reference_start <= region_start and read.reference_end > region_end


class Epipoly:
    def __init__(
        self,
        bam_file,
        bed_file,
        output_tsv,
        check_fasta,
        min_cluster_width,
        max_cluster_width,
        cluster_n,
    ):
        self.bed_path = bed_file
        if check_fasta:
            import pyfaidx

            self.fasta = pyfaidx.Fasta(check_fasta)
        else:
            self.fasta = False
        self.out_tsv = open(output_tsv, "w")
        self.sam = pysam.AlignmentFile(bam_file, "rb")
        self.min_cluster_width = min_cluster_width
        self.max_cluster_width = max_cluster_width
        self.cluster_n = cluster_n

    def iter_bed(self):
        four_cpg = deque(maxlen=self.cluster_n)
        with open(self.bed_path) as bed:
            for line in bed:
                values = line.strip().split("\t")
                try:
                    four_cpg.append((values[0], int(values[1]), int(values[2])))
                except IndexError:
                    raise Exception(
                        "Malformatted bed entry:\n{}"
                        "Each line needs at least three TAB-separated values!".format(
                            line
                        )
                    )
                if len(four_cpg) != self.cluster_n or four_cpg[0][0] != four_cpg[-1][0]:
                    continue
                chrom, start, end = self.update_cpg_locations(four_cpg)
                length = end - start
                if length < self.min_cluster_width or length > self.max_cluster_width:
                    continue
                yield chrom, start, end

    def write_header(self):
        self.out_tsv.write(
            "chr\tstart\tend\tepipolymorphism\tmfrac\tcoverage\tspanning_reads\n"
        )

    def update_cpg_locations(self, cpg_tup):
        self.c_locations = set()
        for cpg in cpg_tup:
            if self.fasta:
                # If the user so desires:
                # check whether the regions in the bed are *really* CpG
                assert self.fasta[cpg[0]][cpg[1]:cpg[2]].seq.upper() == "CG"
            self.c_locations.add(cpg[1])
        clst_chrom = cpg_tup[0][0]
        clst_start = cpg_tup[0][1]
        clst_end = cpg_tup[3][2]
        self.current_bed_entry = clst_chrom, clst_start, clst_end
        return clst_chrom, clst_start, clst_end

    def calc_epi_for_all(self):
        prev_chrom = None
        for chrom, start, stop in self.iter_bed():
            if chrom != prev_chrom:
                print(
                    "Calculating epipolymorphism for CpG sites on "
                    "chromosome {}...".format(chrom)
                )
                prev_chrom = chrom
            reads = self.sam.fetch(contig=chrom, start=start, stop=stop)
            stats_d = self.calc_stats(reads)
            if stats_d["spanning_reads"] > 0:
                out_line = (
                    "{chrom}\t{start}\t{end}\t{epipolymorphism:.5}\t{mfrac:.5}\t"
                    "{coverage}\t{spanning_reads}\n"
                ).format(chrom=chrom, start=start, end=stop, **stats_d)
                self.out_tsv.write(out_line)
                assert stats_d["mfrac"] >= 0 and stats_d["mfrac"] <= 1
                assert (
                    stats_d["epipolymorphism"] >= 0 and stats_d["epipolymorphism"] <= 1
                )

    def calc_stats(self, reads):
        output_dict = {
            "epipolymorphism": "NA",
            "mfrac": "NA",
            "coverage": 0,
            "spanning_reads": 0,
        }
        epi_patterns = []
        n_meth = 0
        for i, read in enumerate(reads):
            # If the read doesn't span the full region, ignore it
            if not spans_region(read, min(self.c_locations), max(self.c_locations) + 1):
                continue

            aln = read.get_aligned_pairs(with_seq=True)
            pos_to_readbase = {x[1]: x[2] for x in aln}

            epi_pattern = []
            for c_loc in self.c_locations:
                # the CpG-dinucleotide as observed in the read (possibly mutated):
                readbase = (pos_to_readbase[c_loc] + pos_to_readbase[c_loc + 1]).upper()
                if readbase not in is_meth:
                    # unexpected dinucleotide (not CG, TG or CA):
                    # there seems to be a SNP, misalignment, or sequencing error,
                    # so let's discard the read
                    epi_pattern = None
                    break
                CpG_is_meth = is_meth[readbase]  # CG = True, TG or CA = False
                epi_pattern.append(CpG_is_meth)
                n_meth += CpG_is_meth

            if epi_pattern:
                epi_patterns.append(tuple(epi_pattern))

        if epi_patterns:
            epi_counts = Counter(epi_patterns)

            # Using the formula of Landan et al. 2012
            n_total = len(epi_patterns)
            epipoly = 1 - sum((n / n_total) ** 2 for n in epi_counts.values())

            # Also calculating the fraction of methylated reads for reference
            mfrac = n_meth / (n_total * self.cluster_n)

            # Returning a dict for easy string formatting
            output_dict["epipolymorphism"] = epipoly
            output_dict["mfrac"] = mfrac
            output_dict["spanning_reads"] = n_total
            output_dict["coverage"] = i + 1
        return output_dict


def main(
    bam_file,
    bed_file,
    output_tsv,
    check_fasta,
    min_cluster_width,
    max_cluster_width,
    cluster_n,
):
    e = Epipoly(
        bam_file,
        bed_file,
        output_tsv,
        check_fasta,
        min_cluster_width,
        max_cluster_width,
        cluster_n,
    )
    e.write_header()
    e.calc_epi_for_all()
    print("Finished writing", output_tsv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculates epipolymorphism as "
        "defined in Landan et al. 2012 from a SAM or BAM file of bisulfite-converted "
        "DNA reads"
    )
    parser.add_argument(
        "bam_file", help="BAM or SAM file of bisulfite-converted, aligned reads"
    )
    parser.add_argument(
        "bed_file",
        help="bed file of genomic CpG sites of interest (usually all CpG sites)",
    )
    parser.add_argument("output_tsv", help="output tsv file path")
    parser.add_argument(
        "--check-fasta",
        help="optionally, use the genome fasta to assert that the CpG locations in "
        "the bed file are really CpG sites",
        default=False,
    )
    parser.add_argument(
        "--min-cluster-width",
        help="CpG clusters below this width will not be reported",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--max-cluster-width",
        help="CpG clusters above this width will not be reported (default: 100; "
        "this should be less or equal to the read length)",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--cluster-n",
        help="number of CpG sites in a cluster (default: 4; this is the usual value)",
        type=int,
        default=4,
    )
    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        main(**vars(args))
