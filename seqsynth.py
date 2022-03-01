#!/usr/bin/env python3
import sys
import os
import gzip
import argparse
import random
import copy
from operator import itemgetter

from utils import reverse_complement, parse_rest, parse_fasta, merge_intervals, get_complement

from colorama import Fore, Back, Style

class SeqSynth:
    def __init__(self, gtf, fa, gene):
        self.gtf  = gtf
        self.fa   = fa
        self.gene = gene
        self.txs = {}
        self.event_spec = []

        # Setup
        self.parse_gtf()
        self.scaffolds = parse_fasta(self.fa)


        '''
        for k, tx in self.txs.items():
            print(k)
            for utr in tx['utrs']:
                print(utr['feature'], utr['strand'])
                utr = self.scaffolds[utr['scaffold']][utr['start']-10:utr['end']+10]
                if 'AATAAA' in utr:
                    print('Found a thing!')
                    print(utr)
                    print(utr.count('AATAAA'))
                if 'TTTATT' in utr:
                    print('Found a reciprocal thing!')
                    print(utr)
                    print(utr.count('TTTATT'))

            print('==========')
        '''

    """
      ____________
     /|           |
    | |  PARSERS  |
    | |___________|
    |/___________/

    """
    def parse_gtf(self):
        print(f'Parsing {Fore.RED}{self.gtf}{Style.RESET_ALL}')
        # Who the hell came up with this stupid file format?
        rest = []
        lb = float('inf')
        ub = 0

        # Check whether annotation is gzipped
        _, ext = os.path.splitext(self.gtf)
        op = lambda x: gzip.open(x, 'rb') if ext == '.gz' else open
        with op(self.gtf) as fh:
            for line in fh:

                # Decode bytestring if .gz
                try:
                    line = line.decode('utf-8')
                except AttributeError as e:
                    pass

                # Skip comments
                if line.startswith('#'):
                    continue
                data = line.split('\t')

                fields = {
                    'scaffold': data[0],
                    'feature': data[2].lower(),
                    'start': int(data[3]) - 1,
                    'end': int(data[4]),
                    'strand': data[6],
                }
                fields['rest'] = parse_rest(data[8], '.gtf')

                if fields['feature'] in ['exon', 'three_prime_utr', 'five_prime_utr']:
                                        #, 'stop_codon', 'start_codon']:
                    gene = fields['rest']['gene_id']
                    tx   = fields['rest']['transcript_id']
                    if gene in self.gene:
                        if tx not in self.txs:
                            self.txs[tx] = {
                                'name': tx,
                                'exons': [],
                                'utrs': [],
                                'scaffold': fields['scaffold'],
                                'strand': fields['strand']
                            }
                        if fields['feature'] == 'exon':
                            fields['rest']['exon_number'] = int(fields['rest']['exon_number'])
                            self.txs[tx]['exons'].append(fields)
                        elif fields['feature'] in ['three_prime_utr', 'utr',
                                                   'five_prime_utr']:
                            self.txs[tx]['utrs'].append(fields)

                        lb = min(fields['start'], lb)
                        ub = max(fields['end'],   ub)

                if fields['feature'] == 'gene' and fields['rest']['gene_id'] in self.gene:
                    lb = min(fields['start'], lb)
                    ub = max(fields['end'],   ub)

        for tx in self.txs.values():
            tx['exons'].sort(key=lambda x: x['rest']['exon_number'])

        # Create introns (maybe move to 'generate_retention')
        for k, tx in self.txs.items():
            ivs = [(f['start'], f['end']) for f in tx['exons']]
            tx['introns'] = get_complement(ivs, lb, ub)

        print(f'Found {len(tx.keys())} transcripts belonging to {Fore.RED}{",".join(self.gene)}{Style.RESET_ALL}')

    def parse_event_spec(self, events):
        types = {'D': 'duplication', 'S': 'skip', 'R': 'retention', 'A': 'APA'}
        for e in events:
            try:
                t = types[e[0].upper()]
                p = float(e[1:])
            except KeyError as err:
                raise KeyError(f'Cannot parse event type {e[0]}')
            except ValueError as err:
                raise ValueError(f'Cannot parse effect size {e[1:]}')

            self.event_spec.append({'type': t, 'pheno': p})

    """
      ________________________
     /|                       |
    | |  SEQUENCE GENERATORS  |
    | |_______________________|
    |/_______________________/

    """
    def generate_rna(self):
        pass

    def generate_duplication(self, tx, dups):
        new_tx = copy.deepcopy(self.txs[tx])
        exons = new_tx['exons']
        i, j = random.choices(range(len(exons)), k=2)
        if j < i:
            i, j = j, i
        exons = exons[:j+1] + exons[i,j+1] + exons[j+1:]
        new_tx['exons'] = exons
        tx_id = f"{new_tx['name']}_D{dups}"
        new_tx['name'] = tx_id
        self.txs[tx_id] = new_tx

    def generate_exon_skip(self, tx, skips):
        new_tx = copy.deepcopy(self.txs[tx])
        exons = new_tx['exons']
        i, j = random.choices(range(len(exons)), k=2)
        if j < i:
            i, j = j, i
        exons = exons[:i] + exons[j+1]
        new_tx['exons'] = exons
        tx_id = f"{new_tx['name']}_S{skips}"
        new_tx['name'] = tx_id
        self.txs[tx_id] = new_tx

    def generate_APA(self, tx, skips):
        # TODO:
        # Figure out APA events

        tx['fwd_apa_motifs'] = sum(x.count('AATAAA') for x in tx['utrs'])
        tx['rv_apa_motifs'] = sum(x.count('TTTATT') for x in tx['utrs'])
        pass

    def generate_retention(self, tx, rets):
        new_tx = copy.deepcopy(self.txs[tx])
        exons = new_tx['exons']
        i = random.choice(range(len(exons)))

        if exons[i]['strand'] == '-':
            pred = lambda x: x[1] <= exons[i]['start']
        else:
            pred = lambda x: x[0] >= exons[i]['end']

        intron = next(x for x in new_tx['introns'] if pred(x))
        intron = {
            'start' : intron[0],
            'end'   : intron[1],
            'strand': exons[i]['strand']
        }

        new_tx['exons'] = exons[:i+1] + [intron] + exons[i+1:]
        tx_id = f"{new_tx['name']}_S{rets}"
        new_tx['name'] = tx_id
        self.txs[tx_id] = new_tx

    def generate_events(self):
        phenotypes = {}
        novel_txs = {}
        dups, skips, apas, rets = 0, 0, 0, 0
        for e in self.event_spec:
            if e['type'] is not None:
                tx = random.choice(self_txs.keys())
                if e['type'] == 'duplication':
                    dups += 1
                    fields = generate_duplication(tx, dups)
                elif e['type'] == 'skip':
                    skips += 1
                    fields = generate_exon_skip(tx, skips)
                elif e['type'] == 'APA':
                    apas += 1
                    fields = generate_APA(tx, apas)
                elif e['type'] == 'retention':
                    rets += 1
                    fields = generate_retention(tx, rets)

                novel_txs[tx] = fields

                if e['pheno'] != 0:
                    phenotypes[tx] = e['pheno']
            else:
                tx = random.choice(self_txs.keys())
                phenotypes[tx] = e['pheno']

        self.txs = dict(self.txs.items(), new_txs.items())
        self.phenotypes = phenotypes

    def find_splice_junctions(self, gtf, fa, od='.', k=31):
        genes = parse_gff(gff_path)
        print('gff parsed')
        scaffolds = parse_fasta(fasta_path)
        print('scaffolds parsed')
        starts = []
        ends = []

        wrong_scaffolds = 0
        for gene, data in genes.items():
            starts = []
            ends = []
            intervals = []
            if data['scaffold'] in scaffolds:
                sequence = scaffolds[data['scaffold']]
            else:
                print(f'{data["scaffold"]} not in FASTA file {fasta_path}')
                wrong_scaffolds += 1
                continue
            for exon in data['exons']:
                intervals.append((int(exon['start'])-1, int(exon['end'])))

            intervals = merge_intervals(intervals)
            sjs = []
            for iv in intervals:
                if iv[1] - iv[0] < 3 * k - 3:
                    # If the length of the exon is < k-1
                    sjs.append(iv)
                else:
                    sjs.append((iv[0], iv[0] + 2 * k - 2))
                    sjs.append((iv[1] - (2 * k) + 2, iv[1]))

            splice_junctions = [
                collapse_N(sequence[iv[0]:iv[1]].upper()) for iv in sjs
            ]

            if data['strand'] == '-':
                splice_junctions = map(reverse_complement, splice_junctions)
            with open(f'{od}/splice_junctions.fasta', 'a') as fh:
                for i, sj in enumerate(splice_junctions):
                    fh.write(f'>{gene}:{i}\n')
                    fh.write(f'{sj}\n')
        print(f'{wrong_scaffolds} scaffolds not found')

    def render_reads(self):
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate simulated RNA reads for a cohort.'
    )
    parser.add_argument('--gtf', help='Transcriptome annotation', type=str)
    parser.add_argument('--fa', help='Genome scaffolds', type=str)
    parser.add_argument('--N', help='Cohort size', type=int)
    parser.add_argument('--depth', help='Mean read depth', type=int)
    parser.add_argument('--genes', help='The gene(s) to generate reads for',
                        nargs='*', type=str)
    parser.add_argument('--events', help='''
            Denote events and phenotypes. D0.1 denotes a duplication event with
            a phenotype 10 percent correlated with it. Other options include A:
            alternative poly-A, S: exon skip, R: intron retention.
            ''', nargs='*', type=str)
    parser.add_argument('--out', help='The output directory', nargs=1, type=str)

    args = parser.parse_args()

    seq = SeqSynth(args.gtf, args.fa, args.genes)
    #seq.parse_event_spec(args.event_spec)

    seq.generate_reads(args.N, args.depth)


