#!/usr/bin/env python3
import sys
import os
import gzip
import argparse
import random
import copy
import uuid
from operator import itemgetter

from utils import reverse_complement, parse_rest, parse_fasta, merge_intervals, get_complement

from numpy import random as nprand
from colorama import Fore, Back, Style

class SeqSynth:
    def __init__(self, gtf, fa, gene):
        self.gtf  = gtf
        self.fa   = fa
        self.gene = gene
        self.txs = {}
        self.novel_txs = {}
        self.event_spec = []
        self.phenotypes = {}

        # Setup
        self.parse_gtf()
        self.scaffolds = parse_fasta(self.fa)

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

        print(f'Found {Fore.GREEN}{len(tx.keys())}{Style.RESET_ALL} transcripts belonging to {Fore.RED}{",".join(self.gene)}{Style.RESET_ALL}')

    def parse_event_spec(self, events):
        types = {
            'D': 'duplication',
            'S': 'skip',
            'R': 'retention',
            'A': 'APA',
            'N': None
        }

        for e in events:
            t, r, p = e.split('/')
            try:
                t = types[t.upper()]
                r = float(r)
                p = [float(p_) for p_ in p.split(':')]
            except KeyError as err:
                raise KeyError(f'Cannot parse event type {t}.')
            except ValueError as err:
                raise ValueError(f'Cannot parse prevalence {r} or effect size {p}.')

            self.event_spec.append({'type': t, 'pheno': p, 'prevalence': r})

    """
      ________________________
     /|                       |
    | |  SEQUENCE GENERATORS  |
    | |_______________________|
    |/_______________________/

    """
    def generate_duplication(self, tx, dups):
        new_tx = copy.deepcopy(self.txs[tx])
        exons = new_tx['exons']
        i, j = random.choices(range(len(exons)), k=2)
        if j < i:
            i, j = j, i
        exons = exons[:j+1] + exons[i:j+1] + exons[j+1:]
        new_tx['exons'] = exons
        tx_id = f"{new_tx['name']}_D{dups}"
        new_tx['name'] = tx_id
        self.novel_txs[tx_id] = new_tx
        print(f'Generated duplication {Fore.RED}{tx_id}{Style.RESET_ALL}')
        return tx_id

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
        self.novel_txs[tx_id] = new_tx
        print(f'Generated exon skip {Fore.RED}{tx_id}{Style.RESET_ALL}')
        return tx_id

    def generate_APA(self, tx, skips):
        # TODO:
        # Figure out APA events

        tx['fwd_apa_motifs'] = sum(x.count('AATAAA') for x in tx['utrs'])
        tx['rv_apa_motifs'] = sum(x.count('TTTATT') for x in tx['utrs'])
        print(f'Generated alternative poly-A {Fore.RED}{tx_id}{Style.RESET_ALL}')
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
            'start'   : intron[0],
            'end'     : intron[1],
            'strand'  : exons[i]['strand'],
            'scaffold': exons[i]['scaffold']
        }

        new_tx['exons'] = exons[:i+1] + [intron] + exons[i+1:]
        tx_id = f"{new_tx['name']}_S{rets}"
        new_tx['name'] = tx_id
        self.novel_txs[tx_id] = new_tx
        print(f'Generated intron retention {Fore.RED}{tx_id}{Style.RESET_ALL}')
        return tx_id

    def generate_events(self):
        dups, skips, apas, rets = 0, 0, 0, 0
        for e in self.event_spec:
            if e['type'] is not None:
                tx = random.choice(list(self.txs.keys()))
                if e['type'] == 'duplication':
                    dups += 1
                    tx_id = self.generate_duplication(tx, dups)
                elif e['type'] == 'skip':
                    skips += 1
                    tx_id = self.generate_exon_skip(tx, skips)
                elif e['type'] == 'APA':
                    apas += 1
                    tx_id = self.generate_APA(tx, apas)
                elif e['type'] == 'retention':
                    rets += 1
                    tx_id = self.generate_retention(tx, rets)

                self.phenotypes[tx_id] = {'effect': e['pheno'],
                                          'prevalence': e['prevalence'],
                                          'type': e['type']}
            else:
                tx = random.choice(self_txs.keys())
                self.phenotypes[tx] = {'effect': e['pheno'],
                                       'prevalence': e['prevalence'],
                                       'type': None}


    """
      _______________________
     /|                      |
    | |  EXPRESSION PROFILE  |
    | |______________________|
    |/______________________/

    """
    def generate_profile(self, N):
        profile = {str(uuid.uuid4())[:7].upper(): [] for _ in range(N)}
        events = copy.deepcopy(self.phenotypes)
        for tx, e in events.items():
            # We want at least one individual with the event, otherwise what is
            # the point?
            e['indiv'] = random.choices(list(profile.keys()),
                                        k=max(1, int(e['prevalence'] * N)))
            for i in e['indiv']:
                profile[i].append(tx)

        self.phenotypes = events
        N_wt_tx = len(self.txs.keys())
        for pn in profile.keys():
            txs = random.choices(list(self.txs.keys()), k=random.choice(range(1, N_wt_tx)))
            profile[pn] = list(set(profile[pn] + txs))
        return profile

    """
      _______________________
     /|                      |
    | |  GENOTYPE RENDERERS  |
    | |______________________|
    |/______________________/

    """
    def render_transcript(self, tx):
        tx = self.txs[tx] if tx in self.txs else self.novel_txs[tx]
        scaffold = self.scaffolds[tx['exons'][0]['scaffold']]
        seq = ''.join(scaffold[e['start']:e['end']] for e in tx['exons'])
        if tx['strand'] == '-':
            seq = reverse_complement(seq)
        return seq

    def render_read(self, tx):
        frag_length = int(nprand.normal(400, 100))
        read_length = int(nprand.normal(100, 25))
        start = max(0, len(tx) - frag_length)
        end = min(len(tx)-1, start+read_length)
        return tx[start:end]

    def render_reads(self, profile, depth, out):
        all_txs = list(self.txs.keys()) + list(self.novel_txs.keys())
        ratios = {tx: random.choice(range(1, 100)) for tx in all_txs}

        txs_seq = {}
        with open(out, 'w') as fh:
            for pn, txs in profile.items():
                denom = sum(ratios[tx] for tx in txs)
                # TODO: implement std dev for expression
                total_expr = nprand.normal(depth, depth/10)
                running_sum = 0

                # R1_fh = open(os.path.join(out, f'{pn}_S1_L001_R1_001.fastq'))
                # R2_fh = open(os.path.join(out, f'{pn}_S1_L001_R2_001.fastq'))
                for tx in txs:
                    if tx not in txs_seq:
                        # If not memoized yet
                        txs_seq[tx] = self.render_transcript(tx)

                    mean = total_expr * ratios[tx]/denom
                    tx_expr = int(nprand.normal(mean, mean/10))
                    reads = (self.render_read(txs_seq[tx]) for _ in range(tx_expr))
                    # for fwd, rev in reads:
                        # render_fastq(pn, fwd, rwd, R1_fh, R2_fh)
                    for i, fwd in enumerate(reads):
                        fh.write(f'>{pn} {tx} {i}')
                        fh.write('\n')
                        fh.write(fwd)
                        fh.write('\n')
                        running_sum += 1

                print(f'Generated {Fore.GREEN}{running_sum}{Style.RESET_ALL} reads of {Fore.GREEN}{len(txs)}{Style.RESET_ALL} transcripts for pn {Fore.RED}{pn}{Style.RESET_ALL} ({Fore.RED}{", ".join(txs)}{Style.RESET_ALL})')

                # close(R1_fh)
                # close(R2_fh)

    """
      _________________________
     /|                        |
    | |  PHENOTYPE GENERATORS  |
    | |________________________|
    |/________________________/

    """
    def generate_phenotypes(self, profile, path):
        try:
            os.makedirs(path)
        except FileExistsError as e:
            pass

        for tx, p in self.phenotypes.items():
            # 1. Generate N samples from a normal distribution.
            # 2. Find pns with correlated phenotype.
            # 3. Multiply the phenotype by the effect size.
            # 4. Write the phenotype to disk.
            pf = p['effect']
            effect = pf[0]
            if len(pf) == 1:
                pf[2] = 1
                pf[1] = effect + 0.5

            while effect <= pf[1]:
                effect_mult = effect + 1.0
                pheno = {k:v for k, v in zip(profile.keys(), nprand.normal(0, 1, len(profile)))}
                for i in p['indiv']:
                    pheno[i] *= effect_mult

                with open(os.path.join(path, f'pheno_{tx}_{effect}'), 'w') as fh:
                    for k,v in pheno.items():
                        fh.write(f'{k}\t{v}\n')

                effect += pf[2]

"""
  ___________
 /|          |
| |  DRIVER  |
| |__________|
|/__________/

"""
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
            Denote events and phenotypes. D/0.04/0.1 denotes a duplication event
            present in 4 percent of the cohort, which correlates 10 percent with
            a phenotype. Other options include A: alternative poly-A, S: exon
            skip, R: intron retention.
            ''', nargs='*', type=str)
    parser.add_argument('--out-fasta', help='The output fasta file', type=str)
    parser.add_argument('--out-pheno', help='The output phenotype file', type=str)

    args = parser.parse_args()

    seq = SeqSynth(args.gtf, args.fa, args.genes)
    seq.parse_event_spec(args.events)
    seq.generate_events()

    profile = seq.generate_profile(args.N)
    seqs = seq.render_reads(profile, args.depth, args.out_fasta)
    seq.generate_phenotypes(profile, args.out_pheno)
