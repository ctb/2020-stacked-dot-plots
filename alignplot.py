#! /usr/bin/env python
"""
Class etc to produce a stacked dotplot and other genome overlap/contamination
stats.
"""
import sys
import argparse
import screed
import matplotlib.pyplot as plt
import csv
import tempfile
import shutil
import subprocess
import os
import glob
from collections import defaultdict


def glob_all(pattern, endings):
    g = []
    for end in endings:
        p = pattern + end
        g.extend(glob.glob(p))

    return g


class StackedDotPlot:
    """
    Build a stacked dot plot.

    Takes:
    * query accession,
    * multiple target accessions,
    * an optional info file containing mappings from accession to names
    * an optional directory containing the genomes
    """
    endings = '.gz', '.fa', '.fna'
    use_mashmap = False

    def __init__(self, q_acc, t_acc_list, info_file=None, genomes_dir=None):
        if genomes_dir is None:
            genomes_dir = '.'
        else:
            genomes_dir = genomes_dir.rstrip('/')


        queryfiles = glob_all(f'{genomes_dir}/{q_acc}*', self.endings)
        print(queryfiles)
        assert len(queryfiles) == 1, queryfiles
        queryfile = queryfiles[0]
        print(f'found queryfile for {q_acc}: {queryfile}')
        self.queryfile = queryfile
        self.q_acc = q_acc

        targetfiles = []
        for t_acc in t_acc_list:
            x = glob_all(f'{genomes_dir}/{t_acc}*', self.endings)
            assert len(x) == 1, x
            targetfiles.append(x[0])
            print(f'found targetfile for {t_acc}: {x[0]}')

        self.targetfiles = targetfiles
        self.t_acc_list = list(t_acc_list)
        
        self.query_name = q_acc
        self.target_names = {}
        for acc in t_acc_list:
            self.target_names[acc] = acc

        if info_file:
            for row in csv.DictReader(open(info_file, 'rt')):
                if self.q_acc == row['acc']:
                    self.query_name = row['ncbi_tax_name']
                if row['acc'] in self.target_names:
                    self.target_names[row['acc']] = row['ncbi_tax_name']

        self.q_starts = {}
        self.q_sofar = 0

    def __call__(self,):
        "Run all the things, produce a plot."
        results = {}

        for t_acc, targetfile in zip(self.t_acc_list, self.targetfiles):
            name = self.target_names[t_acc]

            if self.use_mashmap:
                sum_shared, starts, ends = self.run_mashmap(targetfile)
            else:
                sum_shared, starts, ends = self.run_nucmer(targetfile)

            print(f'shared w/{name}: {sum_shared:.1f}kb')
            results[t_acc] = (sum_shared, starts, ends)

        self.results = results
        return self.plot()

    def run_mashmap(self, targetfile):
        "Run mashmap. Deprecated."
        print("running mashmap...")
        tempdir = tempfile.mkdtemp()
        outfile = os.path.join(tempdir, "mashmap.out")
        cmd = f"mashmap -q {self.queryfile} -r {targetfile} -o {outfile} --pi 95" # -f none -s 1000
        print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)
        
        print(f"...done! reading output from {outfile}.")

        results = self._read_mashmap(outfile)
        shutil.rmtree(tempdir)
        return results
    
    def _read_mashmap(self, filename):
        "Parse the mashmap output."
        fp = open(filename, 'rt')
        starts = []
        ends = []

        q_starts = self.q_starts
        q_sofar = self.q_sofar
        t_starts = {}
        t_sofar = 0
        
        sumshared = 0
        matches = defaultdict(list)
        for line in fp:
            line = line.strip().split()
            query, qsize, qstart, qend, strand, target, tsize, tstart, tend, pident = line
            qsize = int(qsize) / 1e3
            qstart = int(qstart) / 1e3
            qend = int(qend) / 1e3
            tsize = int(tsize) / 1e3
            tstart = int(tstart) / 1e3
            tend = int(tend) / 1e3
            pident = float(pident)

            q_base = q_starts.get(query)
            if q_base is None:
                size = qsize
                q_starts[query] = q_sofar
                q_base = q_sofar
                q_sofar += size

            t_base = t_starts.get(target)
            if t_base is None:
                size = tsize
                t_starts[target] = t_sofar
                t_base = t_sofar
                t_sofar += size

            assert qend > qstart
            sumshared += qend - qstart

            starts.append((t_base + tstart, q_base + qstart))
            ends.append((t_base + tend, q_base + qend))

        self.q_sofar = q_sofar
        return sumshared, starts, ends
        
    def run_nucmer(self, targetfile):
        "Run nucmer and show coords."
        print("running nucmer & show-coords...")
        tempdir = tempfile.mkdtemp()

        queryfile = self.queryfile
        if self.queryfile.endswith('.gz'):
            queryfile = os.path.join(tempdir, "query.fa")
            subprocess.check_call(f"gunzip -c {self.queryfile} > {queryfile}", shell=True)

        if targetfile.endswith('.gz'):
            newfile = os.path.join(tempdir, "target.fa")
            subprocess.check_call(f"gunzip -c {targetfile} > {newfile}", shell=True)
            targetfile = newfile
            
        cmd = f"nucmer -p {tempdir}/cmp {queryfile} {targetfile} 2> /dev/null"
        #print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        deltafile = f"{tempdir}/cmp.delta"
        coordsfile = f"{tempdir}/cmp.coords"

        cmd = f"show-coords -T {deltafile} > {coordsfile} 2> /dev/null"
        #print(f"running {cmd}")
        subprocess.check_call(cmd, shell=True)

        print(f"...done! reading output from {tempdir}.")

        results = self._read_nucmer(coordsfile)
        #shutil.rmtree(tempdir)
        return results
    
    def _read_nucmer(self, filename):
        "Parse the nucmer output."
        fp = open(filename, 'rt')
        lines = fp.readlines()
        assert lines[1].startswith('NUCMER'), (filename, lines[0])
        assert not lines[2].strip()

        sumshared = 0
        starts = []
        ends = []
        q_starts = self.q_starts
        q_sofar = self.q_sofar
        t_starts = {}
        t_sofar = 0
        
        for line in lines[4:]:
            line = line.strip().split('\t')
            qstart, qend, tstart, tend, qsize, tsize, pident, query, target = line
            qsize = int(qsize) / 1e3
            qstart = int(qstart) / 1e3
            qend = int(qend) / 1e3
            tsize = int(tsize) / 1e3
            tstart = int(tstart) / 1e3
            tend = int(tend) / 1e3
            pident = float(pident)
  
            # identity and length filter
            if pident < 95 or abs(qend - qstart) < 0.5:
                continue

            q_base = q_starts.get(query)
            if q_base is None:
                size = qsize
                q_starts[query] = q_sofar
                q_base = q_sofar
                q_sofar += size

            t_base = t_starts.get(target)
            if t_base is None:
                size = tsize
                t_starts[target] = t_sofar
                t_base = t_sofar
                t_sofar += size

            sumshared += int(qend) - int(qstart)

            starts.append((t_base + tstart, q_base + qstart))
            ends.append((t_base + tend, q_base + qend))

        self.q_sofar = q_sofar
        return sumshared, starts, ends
        
    def plot(self):
        "Do the actual stacked dotplot plotting."
        plt.ylabel(f'{self.q_acc}: {self.query_name}')

        colors = ('r-', 'b-', 'g-')

        # the use of max_x is what makes it a stacked dotplot!! :)
        max_x = 0                         # track where to start each target

        # iterate over each set of features, plotting lines.
        for t_acc, color in zip(self.t_acc_list, colors):
            name = self.target_names[t_acc]
            sum_shared, starts, ends = self.results[t_acc]

            this_max_x = 0
            line = None
            for (x_0, y_0), (x_1, y_1) in zip(starts, ends):
                line = plt.plot((x_0 + max_x, x_1 + max_x), (y_0, y_1), color)
                this_max_x = max(this_max_x, x_0, x_1)

            if line:
                line[0].set_label(name)

            # "stack" the dotplots horizontally.
            max_x = this_max_x

        plt.legend(loc='lower right')
            
        return plt.gcf()


def main():
    dotplot = StackedDotPlot('GCA_003222275.1',
                             ['GCA_003220225.1', 'GCA_003222275.1'],
                             'list.csv', './genomes')
    _ = dotplot()

    print('saving')
    plt.savefig('/tmp/test.png')

    dotplot.use_mashmap = True
    _ = dotplot()

    print('saving')
    plt.savefig('/tmp/test-mashmap.png')
    return 0


if __name__ == '__main__':
    sys.exit(main())
