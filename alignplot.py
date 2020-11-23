#! /usr/bin/env python
"""
Class etc to produce a stacked dotplot and other genome overlap/contamination
stats.

TODO:
* argparse the thang
"""
import sys
import argparse
import matplotlib.pyplot as plt
import csv
import tempfile
import shutil
import subprocess
import os
import glob
from collections import defaultdict, namedtuple
import numpy

import screed
from interval import interval


AlignedRegion = namedtuple(
    "AlignedRegion", "query, target, qstart, qend, tstart, tend, pident, qsize, tsize"
)


def glob_all(pattern, endings):
    g = []
    for end in endings:
        p = pattern + end
        g.extend(glob.glob(p))

    return g


def group_regions_by(regions, by_type):
    "gather regions by 'target' or 'query'"

    assert by_type in ("target", "query")

    regions_by = defaultdict(list)
    for region in regions:
        name = getattr(region, by_type)
        x = regions_by[name]
        x.append(region)

    return regions_by


def calc_regions_aligned_bp(regions_by, by_type, filter_by=None):
    assert by_type in ("target", "query")

    regions_sum_kb = {}
    for name, rr in regions_by.items():
        if filter_by:
            rr = [r for r in rr if filter_by(r)]
        if not rr:
            continue

        if by_type == "target":
            ivals = [interval[r.tstart, r.tend] for r in rr]
        elif by_type == "query":
            ivals = [interval[r.qstart, r.qend] for r in rr]
        else:
            raise Exception(by_type)

        ii = interval()
        for iv in ivals:
            ii |= iv

        sum_kb = 0
        for component in ii:
            (start, end) = component
            size = end - start
            sum_kb += size

        regions_sum_kb[name] = sum_kb
    return regions_sum_kb


def region_size(region, by_type):
    assert by_type in ("target", "query")

    if by_type == "target":
        return abs(region.tend - region.tstart)
    elif by_type == "query":
        return abs(region.qend - region.qstart)
    raise Exception(f"unhandled by_type {by_type}")


def load_contig_sizes(genomefile):
    "load in all the actual contig sizes for this genome (in kb)"
    all_sizes = {}
    for record in screed.open(genomefile):  # @CTB
        all_sizes[record.name.split()[0]] = len(record.sequence) / 1e3

    return all_sizes


class StackedDotPlot:
    """
    Build a stacked dot plot.

    Takes:
    * query accession,
    * multiple target accessions,
    * an optional info file containing mappings from accession to names
    * an optional directory containing the genomes
    """

    endings = ".gz", ".fa", ".fna"

    def __init__(
        self, q_acc, t_acc_list, info_file=None, genomes_dir=None, use_mashmap=False
    ):
        self.use_mashmap = use_mashmap

        if genomes_dir is None:
            genomes_dir = "."
        else:
            genomes_dir = genomes_dir.rstrip("/")

        queryfiles = glob_all(f"{genomes_dir}/{q_acc}*", self.endings)
        print(queryfiles)
        assert len(queryfiles) == 1, queryfiles
        queryfile = queryfiles[0]
        print(f"found queryfile for {q_acc}: {queryfile}")
        self.queryfile = queryfile
        self.q_acc = q_acc

        targetfiles = []
        for t_acc in t_acc_list:
            x = glob_all(f"{genomes_dir}/{t_acc}*", self.endings)
            assert len(x) == 1, x
            targetfiles.append(x[0])
            print(f"found targetfile for {t_acc}: {x[0]}")

        self.targetfiles = targetfiles
        self.t_acc_list = list(t_acc_list)

        self.query_name = q_acc
        self.target_names = {}
        for acc in t_acc_list:
            self.target_names[acc] = acc

        if info_file:
            for row in csv.DictReader(open(info_file, "rt")):
                if self.q_acc == row["acc"]:
                    self.query_name = row["ncbi_tax_name"]
                if row["acc"] in self.target_names:
                    self.target_names[row["acc"]] = row["ncbi_tax_name"]

    def get_targetfile(self, t_acc):
        targetfile = None
        for find_t_acc, find_targetfile in zip(self.t_acc_list, self.targetfiles):
            if find_t_acc == t_acc:
                targetfile = find_targetfile
                break

        assert targetfile
        return targetfile

    def __call__(self):
        "Run all the things, produce a plot."
        results = {}

        for t_acc, targetfile in zip(self.t_acc_list, self.targetfiles):
            name = self.target_names[t_acc]

            if self.use_mashmap:
                regions = self.run_mashmap(targetfile)
            else:
                regions = self.run_nucmer(targetfile)

            results[t_acc] = regions

        self.results = results
        return self.plot()

    def run_mashmap(self, targetfile):
        "Run mashmap. Deprecated."
        print("running mashmap...")
        tempdir = tempfile.mkdtemp()
        outfile = os.path.join(tempdir, "mashmap.out")
        cmd = f"mashmap -q {self.queryfile} -r {targetfile} -o {outfile} --pi 95"  # -f none -s 1000
        print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        print(f"...done! reading output from {outfile}.")

        results = self._read_mashmap(outfile)
        shutil.rmtree(tempdir)
        return results

    def _read_mashmap(self, filename):
        "Parse the mashmap output."
        fp = open(filename, "rt")

        regions = []
        for line in fp:
            line = line.strip().split()
            (
                query,
                qsize,
                qstart,
                qend,
                strand,
                target,
                tsize,
                tstart,
                tend,
                pident,
            ) = line
            region = AlignedRegion(
                qsize=int(qsize) / 1e3,
                qstart=int(qstart) / 1e3,
                qend=int(qend) / 1e3,
                tsize=int(tsize) / 1e3,
                tstart=int(tstart) / 1e3,
                tend=int(tend) / 1e3,
                pident=float(pident),
                query=query,
                target=target,
            )

            assert region.qend > region.qstart
            regions.append(region)

        return regions

    def run_nucmer(self, targetfile):
        "Run nucmer and show coords."
        print(f"running nucmer & show-coords for {targetfile}...")
        tempdir = tempfile.mkdtemp()

        queryfile = self.queryfile
        if self.queryfile.endswith(".gz"):
            queryfile = os.path.join(tempdir, "query.fa")
            subprocess.check_call(
                f"gunzip -c {self.queryfile} > {queryfile}", shell=True
            )

        if targetfile.endswith(".gz"):
            newfile = os.path.join(tempdir, "target.fa")
            subprocess.check_call(f"gunzip -c {targetfile} > {newfile}", shell=True)
            targetfile = newfile

        cmd = f"nucmer -p {tempdir}/cmp {queryfile} {targetfile} 2> /dev/null"
        # print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        deltafile = f"{tempdir}/cmp.delta"
        coordsfile = f"{tempdir}/cmp.coords"

        cmd = f"show-coords -T {deltafile} > {coordsfile} 2> /dev/null"
        # print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        print(f"...done! reading output from {tempdir}.")

        results = self._read_nucmer(coordsfile)
        # shutil.rmtree(tempdir)
        return results

    def _read_nucmer(self, filename):
        "Parse the nucmer output."
        fp = open(filename, "rt")
        lines = fp.readlines()
        assert lines[1].startswith("NUCMER"), (filename, lines[0])
        assert not lines[2].strip()

        regions = []
        for line in lines[4:]:
            line = line.strip().split("\t")
            qstart, qend, tstart, tend, qsize, tsize, pident, query, target = line
            region = AlignedRegion(
                qsize=int(qsize) / 1e3,
                qstart=int(qstart) / 1e3,
                qend=int(qend) / 1e3,
                tsize=int(tsize) / 1e3,
                tstart=int(tstart) / 1e3,
                tend=int(tend) / 1e3,
                pident=float(pident),
                query=query,
                target=target,
            )

            # @CTB check orientation somehow!

            # identity and length filter - @CTB move outside!
            #            if region.pident < 95 or abs(region.qend - region.qstart) < 0.5:
            #                continue

            regions.append(region)

        return regions

    def plot(self):
        "Do the actual stacked dotplot plotting."
        if self.q_acc == self.query_name:
            ylabel_text = f"{self.q_acc} (coords in kb)"
        else:
            ylabel_text = f"{self.q_acc}: {self.query_name} (kb)"
        plt.ylabel(ylabel_text)
        plt.xlabel("coordinates of matches (scaled to kb)")

        colors = ("r-", "b-", "g-")

        q_starts = {}
        q_sofar = 0

        # the use of max_x is what makes it a stacked dotplot!! :)
        max_x = 0  # track where to start each target

        # iterate over each set of features, plotting lines.
        for t_acc, color in zip(self.t_acc_list, colors):
            name = self.target_names[t_acc]
            # @CTB if we move this out of the loop and plot self-x-self
            # there is an interestng effect of showing distribution. exploreme!
            t_starts = {}
            t_sofar = 0

            sum_shared = 0
            line = None
            this_max_x = 0
            for region in self.results[t_acc]:
                sum_shared += region.qend - region.qstart

                # calculate the base y position for this query contig --
                q_base = q_starts.get(region.query)
                if q_base is None:
                    q_starts[region.query] = q_sofar
                    q_base = q_sofar
                    q_sofar += region.qsize

                # calculate the base x position for this target contig --
                t_base = t_starts.get(region.target)
                if t_base is None:
                    t_starts[region.target] = t_sofar
                    t_base = t_sofar
                    t_sofar += region.tsize

                x_0 = t_base + region.tstart
                y_0 = q_base + region.qstart

                x_1 = t_base + region.tend
                y_1 = q_base + region.qend

                # stack 'em horizontally with max_x
                line = plt.plot((x_0 + max_x, x_1 + max_x), (y_0, y_1), color)
                this_max_x = max(this_max_x, x_0, x_1)

            # label the last plotted line w/the right name to make legend
            if line:
                line[0].set_label(name)

            # "stack" the dotplots horizontally.
            max_x = this_max_x
            print(f"shared w/{name}: {sum_shared:.1f}kb")

        plt.legend(loc="lower right")

        return plt.gcf()

    def target_response_curve(self, t_acc):
        regions = self.results[t_acc]

        # first, find the targetfile (genome) for this accession
        targetfile = self.get_targetfile(t_acc)

        # calculate and sort region summed kb in alignments over 95%
        regions_by_target = group_regions_by(regions, "target")
        regions_aligned_kb = calc_regions_aligned_bp(
            regions_by_target, "target", filter_by=lambda r: r.pident >= 95
        )
        region_items = list(regions_aligned_kb.items())
        region_items.sort(key=lambda x: -x[1])

        # load in all the actual contig sizes for this genome
        all_sizes = load_contig_sizes(targetfile)
        sum_bp = sum(all_sizes.values())

        # construct points for plot --
        x = [0]  # kb in target contigs
        y = [0]  # alignments in those target contigs
        sofar = 0
        aligned_sofar = 0

        # start with contigs with most aligned bases first - the sorting order matters here!
        for name, ani95_kb in region_items:
            # ok, track total kb and aligned kb added by this contig
            sofar += all_sizes[name]
            aligned_sofar += ani95_kb
            assert all_sizes[name] > 0

            x.append(sofar)
            y.append(aligned_sofar)

        saturation_point = sofar

        # add in the rest of the contigs that have no alignments in 'em'
        remaining_names = set(all_sizes) - set(region_items)
        for contig in remaining_names:
            sofar += all_sizes[contig]
            x.append(sofar)
            y.append(aligned_sofar)

        return numpy.array(x), numpy.array(y), saturation_point

    def query_response_curve(self):
        # aggregate regions over _all_ results
        regions = []
        for k, v in self.results.items():
            regions.extend(v)

        queryfile = self.queryfile

        # calculate and sort region summed kb in alignments over 95%
        regions_by_query = group_regions_by(regions, "query")
        regions_aligned_kb = calc_regions_aligned_bp(
            regions_by_query, "query", filter_by=lambda r: r.pident >= 95
        )
        region_items = list(regions_aligned_kb.items())
        region_items.sort(key=lambda x: -x[1])

        # load in all the actual contig sizes for this genome
        all_sizes = load_contig_sizes(queryfile)
        sum_bp = sum(all_sizes.values())

        # construct points for plot --
        x = [0]  # kb in query contigs
        y = [0]  # alignments in those query contigs
        sofar = 0
        aligned_sofar = 0

        # start with contigs with most aligned bases first - the sorting order matters here!
        for name, ani95_kb in region_items:
            # ok, track total kb and aligned kb added by this contig
            sofar += all_sizes[name]
            aligned_sofar += ani95_kb
            assert all_sizes[name] > 0

            x.append(sofar)
            y.append(aligned_sofar)

        saturation_point = sofar

        # add in the rest of the contigs that have no alignments in 'em'
        remaining_names = set(all_sizes) - set(region_items)
        for contig in remaining_names:
            sofar += all_sizes[contig]
            x.append(sofar)
            y.append(aligned_sofar)

        return numpy.array(x), numpy.array(y), saturation_point


def main():
    p = argparse.ArgumentParser()
    p.add_argument("query_acc")
    p.add_argument("target_accs", nargs="+")
    p.add_argument(
        "-g",
        "--genomes-directory",
        default="./genomes",
        help="directory with genome files in it",
    )
    p.add_argument("-i", "--info-file", help="CSV file with nicer names for accessions")
    p.add_argument("-o", "--output-prefix", default="alignplot-")
    args = p.parse_args()

    dotplot = StackedDotPlot(
        args.query_acc, args.target_accs, args.info_file, args.genomes_directory
    )
    _ = dotplot()

    print(f"saving {args.output_prefix}-nucmer.png")
    plt.savefig(f"{args.output_prefix}-nucmer.png")
    plt.cla()

    dotplot.use_mashmap = True
    _ = dotplot()

    print(f"saving {args.output_prefix}-mashmap.png")
    plt.savefig(f"{args.output_prefix}-mashmap.png")
    plt.cla()

    t_acc = dotplot.t_acc_list[0]
    x, y, sat1 = dotplot.target_response_curve(t_acc)
    x2, y2, sat2 = dotplot.query_response_curve()

    plt.plot(x, y / max(y), "b-", label=f"target loss ({t_acc})")
    plt.plot(x2, y2 / max(y2), "g-", label="query loss")

    plt.xlabel("kb in genome contigs removed")
    plt.ylabel("fraction of alignments removed")
    plt.legend(loc="lower right")

    print(f"saving {args.output_prefix}-response.png")
    plt.savefig(f"{args.output_prefix}-response.png")
    plt.cla()

    return 0


if __name__ == "__main__":
    sys.exit(main())
