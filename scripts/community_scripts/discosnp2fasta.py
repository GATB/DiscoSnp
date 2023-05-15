#! /usr/bin/env python3

#######################################
# A script to generated fasta, nexus and phylip files from 
# DiscoRad/DiscoRadSnp generated pseudo-fasta files
# Also generates multiline fasta like file (.pseudo.fasta)
# Also generates partition list
#
# depends on Biopython and Pandas
# 
# The script filters by max missing data per genotype (missingness), minimum polymorphic sites, and rank. It uses a lookup table to change the standard sample names (G1, G2, etc) to actual sample names. It has the option to include locus name in the header (defline) of the pseudo-fasta file (every block of entries is a different locus). It outputs fasta, nexus and phylip files either as alleles, or as genotypes (collapsed alleles).
#
# The script  can be run, for example, as:
# discosnp2fasta.py -i discoRad_k_31_c_3_D_30_P_5_m_5_raw_filtered.fa -o callitrichids -l sample_correspondence -max_miss .1 -min_poly 5
#
# Only the -i (input file), -o (base name of output file) and -l (tab separated G to sample look up table) are required. All other parameters have reasonable defaults, including -path (which assumes current directory).
#
# Author:  Tomas Hrbek
# Email: hrbek@evoamazon.net
# Date: 20.04.2023
# Version: 0.9
#######################################

__author__ = 'legal'

import os
import re
import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Script to generate a fasta file from DiscoSnpRad/DiscoSnp++ pseudo-fasta file')
parser.add_argument('-i','--infile', help='DiscoSnpRad/DiscoSnp++ pseudo-fasta file', type=str, required=True)
parser.add_argument('-o','--outfile', help='Output file base name (fasta and nexus)', type=str, required=True)
parser.add_argument('-l','--lookup', help='Gx to sample lookup table (tab separated); in project directory', type=str, required=True)
parser.add_argument('-path', help='Path to project directory; default ./', type=str, required=False, default='./')
parser.add_argument('-cons', help='Generate consensus sequence; default True', type=str, required=False, default='True')
parser.add_argument('-loc_info', help='Include locus name in header; default False', type=str, required=False, default='False')
parser.add_argument('-min_rank','--minimum_rank', help='minimum rank; parallog metric, default .4 (decimal)', type=float, required=False, default=.4)
parser.add_argument('-max_miss','--maximum_missingness', help='maximum missing data per locus, default .5 (decimal)', type=float, required=False, default=.5)
parser.add_argument('-min_poly','--minimum_polymorphisms', help='minimum number of SNPs per locus, default 3 (integer)', type=int, required=False, default=3)
args = parser.parse_args()

fasta_in = args.infile
fasta_out = args.outfile
lookup = args.lookup
path = args.path
# solving issues of parsing boolians: https://from-locals.com/python-argparse-bool/
if args.cons == 'True' or args.cons == 'TRUE' or args.cons == 'T' :
    cons = True
else:
    cons = False
if args.loc_info == 'True' or args.loc_info == 'TRUE' or args.loc_info == 'T' :
    loc = True
else:
    loc = False
min_rank = args.minimum_rank
max_miss = args.maximum_missingness
min_poly = args.minimum_polymorphisms


# generate id reads keys
columns = ['id','taxon']
df = pd.read_csv(os.path.join(path, lookup), header=None, names=columns, sep='\t')
r = dict(zip(df['id'], df['taxon']))
# generate a key for holding concatenated sequences (s = consensus; s0&s1 = alleles)
s = [' '] * len(df['id'])
s = dict(zip(df['taxon'], s))
s0 = [' '] * len(df['id'])
s0 = dict(zip(df['taxon'], s0))
s1 = [' '] * len(df['id'])
s1 = dict(zip(df['taxon'], s1))
# generate a key holding name lengths
taxon_len = [len(x) for x in df['taxon']]
taxon_pad = [max(taxon_len) + 5 - x for x in taxon_len]
pad = dict(zip(df['taxon'], taxon_pad))


def extract_fasta(seq, sample, allele, loc):
    locus_pattern = re.compile(r"(.+?)(?=\|)")
    # extract locus name
    locus = locus_pattern.match(seq.id).group(1)
    # change fasta header ID
    if loc == True:
        seq_id = '>' + r.get(sample.group(1), 'NA') + '_' + allele + '_' + locus
    else:
        seq_id = '>' + r.get(sample.group(1), 'NA') + '_' + allele
    if sample.group(2) != '0' and sample.group(2) != '1':
        seq_seq = ''.join('N' for i in range(len(seq.seq)))
    else:
        seq_seq = seq.seq
    return sample.group(1), seq_id, seq_seq;
    
def consense_fasta(seq_a, seq_b, sample, allele, loc):
    locus_pattern = re.compile(r"(.+?)(?=\|)")
    # extract locus name
    locus = locus_pattern.match(seq_a.id).group(1)
    # change fasta header ID
    if loc == True:
        seq_id = '>' + r.get(sample.group(1), 'NA') + '_' + allele + '_' + locus
    else:
        seq_id = '>' + r.get(sample.group(1), 'NA') + '_' + locus
    if sample.group(2) == '0' and sample.group(3) == '0':
        seq_seq = seq_a.seq
    elif sample.group(2) == '1' and sample.group(3) == '1':
        seq_seq = seq_b.seq
    elif sample.group(2) != '0' and sample.group(2) != '1':
        seq_seq = ''.join('N' for i in range(len(seq_a.seq)))
    else:
        seq = []
        for pos in range (0, min(len(seq_a.seq), len(seq_b.seq))):
            if seq_a.seq[pos] == seq_b.seq[pos]:
                seq.append(seq_a.seq[pos])
            elif (seq_a.seq[pos] == 'G' and seq_b.seq[pos] == 'A') or (seq_a.seq[pos] == 'A' and seq_b.seq[pos] == 'G'):
                seq.append('R')
            elif (seq_a.seq[pos] == 'C' and seq_b.seq[pos] == 'T') or (seq_a.seq[pos] == 'T' and seq_b.seq[pos] == 'C'):
                seq.append('Y')
            elif (seq_a.seq[pos] == 'G' and seq_b.seq[pos] == 'C') or (seq_a.seq[pos] == 'C' and seq_b.seq[pos] == 'G'):
                seq.append('S')
            elif (seq_a.seq[pos] == 'A' and seq_b.seq[pos] == 'T') or (seq_a.seq[pos] == 'T' and seq_b.seq[pos] == 'A'):
                seq.append('W')
            elif (seq_a.seq[pos] == 'G' and seq_b.seq[pos] == 'T') or (seq_a.seq[pos] == 'T' and seq_b.seq[pos] == 'G'):
                seq.append('K')
            elif (seq_a.seq[pos] == 'C' and seq_b.seq[pos] == 'A') or (seq_a.seq[pos] == 'A' and seq_b.seq[pos] == 'C'):
                seq.append('M')
            # if indel heterozygote, makes the genotype a transversion heterozygote
            # other option make indel heterozygote a homozygote
            # ideally indels will have been filtered out a priori
            elif (seq_a.seq[pos] == '-' and seq_b.seq[pos] == 'A') or (seq_a.seq[pos] == 'A' and seq_b.seq[pos] == '-'):
                seq.append('W')
            elif (seq_a.seq[pos] == '-' and seq_b.seq[pos] == 'G') or (seq_a.seq[pos] == 'G' and seq_b.seq[pos] == '-'):
                seq.append('S')
            elif (seq_a.seq[pos] == '-' and seq_b.seq[pos] == 'C') or (seq_a.seq[pos] == 'C' and seq_b.seq[pos] == '-'):
                seq.append('M')
            elif (seq_a.seq[pos] == '-' and seq_b.seq[pos] == 'T') or (seq_a.seq[pos] == 'T' and seq_b.seq[pos] == '-'):
                seq.append('K')
            else:
                seq.append('N')
        seq_seq = ''.join(str(x) for x in seq)
    return sample.group(1), seq_id, seq_seq;

def make_fasta(s, path, fasta_out):
    with open(os.path.join(path, fasta_out + '.fas'), 'w') as output_handle:
        for x in s.keys():
            output_handle.write(''.join('>' + x + '\n'))
            output_handle.write(''.join(s.get(x) + '\n'))

def make_fasta_alleles(s0, s1, path, fasta_out):
    with open(os.path.join(path, fasta_out + '.fas'), 'w') as output_handle:
        for x in s0.keys():
            output_handle.write(''.join('>' + x + '_0\n'))
            output_handle.write(''.join(s0.get(x) + '\n'))
            output_handle.write(''.join('>' + x + '_1\n'))
            output_handle.write(''.join(s1.get(x) + '\n'))

def make_nexus(s, path, fasta_out):
    with open(os.path.join(path, fasta_out + '.nex'), 'w') as output_handle:
        output_handle.write(''.join('#NEXUS\n'))
        output_handle.write(''.join('BEGIN DATA;\n'))
        output_handle.write(''.join('  DIMENSIONS NTAX={0} NCHAR={1};\n'.format(len(s.keys()), len(s.get(list(s.keys())[0])))))
        output_handle.write(''.join('  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=NO;\n'))
        output_handle.write(''.join('  MATRIX\n'))
        for x in s.keys():
            output_handle.write(''.join(x))
            output_handle.write(''.join(' ' for i in range(pad.get(x))))
            output_handle.write(''.join(s.get(x) + '\n'))
        output_handle.write(''.join(';\n'))
        output_handle.write(''.join('END;\n\n'))
        # add partition information
        output_handle.write(''.join('BEGIN SETS;\n'))
        with open(os.path.join(path, fasta_out + '.partitions'), 'r') as partition_handle:
            output_handle.write(''.join(partition_handle.read()))
        output_handle.write(''.join('END;\n'))

def make_nexus_alleles(s0, s1, path, fasta_out):
    with open(os.path.join(path, fasta_out + '.nex'), 'w') as output_handle:
        output_handle.write(''.join('#NEXUS\n'))
        output_handle.write(''.join('BEGIN DATA;\n'))
        output_handle.write(''.join('  DIMENSIONS NTAX={0} NCHAR={1};\n'.format(len(s0.keys()), len(s0.get(list(s0.keys())[0])))))
        output_handle.write(''.join('  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=NO;\n'))
        output_handle.write(''.join('  MATRIX\n'))
        for x in s0.keys():
            output_handle.write(''.join(x + '_0'))
            output_handle.write(''.join(' ' for i in range(pad.get(x))))
            output_handle.write(''.join(s0.get(x) + '\n'))
            output_handle.write(''.join(x + '_1'))
            output_handle.write(''.join(' ' for i in range(pad.get(x))))
            output_handle.write(''.join(s1.get(x) + '\n'))
        output_handle.write(''.join(';\n'))
        output_handle.write(''.join('END;\n\n'))
        # add partition information
        output_handle.write(''.join('BEGIN SETS;\n'))
        with open(os.path.join(path, fasta_out + '.partitions'), 'r') as partition_handle:
            output_handle.write(''.join(partition_handle.read()))
        output_handle.write(''.join('END;\n'))

def make_phylip(s, path, fasta_out):
    with open(os.path.join(path, fasta_out + '.phy'), 'w') as output_handle:
        output_handle.write(''.join('{0} {1}\n'.format(str(int(len(s.keys()))), len(s.get(list(s.keys())[0])))))
        for x in s.keys():
            output_handle.write(''.join(x))
            output_handle.write(''.join(' ' for i in range(pad.get(x))))
            output_handle.write(''.join(s.get(x) + '\n'))

def make_phylip_alleles(s0, s1, path, fasta_out):
    with open(os.path.join(path, fasta_out + '.phy'), 'w') as output_handle:
        output_handle.write(''.join('{0} {1}\n'.format(str(int(len(s0.keys())*2)), len(s0.get(list(s0.keys())[0])))))
        for x in s0.keys():
            output_handle.write(''.join(x + '_0'))
            output_handle.write(''.join(' ' for i in range(pad.get(x))))
            output_handle.write(''.join(s0.get(x) + '\n'))
            output_handle.write(''.join(x + '_1'))
            output_handle.write(''.join(' ' for i in range(pad.get(x))))
            output_handle.write(''.join(s1.get(x) + '\n'))

########
# main
with open(os.path.join(path, fasta_out + '.pseudo.fasta'), 'w') as output_handle:
    with open(os.path.join(path, fasta_out + '.partitions'), 'w') as partition_handle:
        with open(os.path.join(path, fasta_in), 'r') as input_handle:
            records = list(SeqIO.parse(input_handle, 'fasta'))
            i = 0
            j = 0
            l = 0
            sample_pattern = re.compile(r"(G[0-9]+)_([01\.])/([01\.])")
            missing_pattern = re.compile(r"G[0-9]+_\./\.")
            #poly_pattern = re.compile(r"P_[0-9]+:[0-9]+_[CTAG]/[CTAG]")
            poly_pattern = re.compile(r".+nb_pol_([0-9]+)")
            rank_pattern = re.compile(r".+rank_([\.0-9]+)")
            while i < len(records):
                #filter by rank
                rank = float(rank_pattern.match(records[i].id).group(1))
                if rank < min_rank:
                    i += 2
                    continue
                # filter on number of polymorphic sites per locus
                poly = int(poly_pattern.match(records[i].id).group(1))
                if poly < min_poly:
                    i += 2
                    continue
                # filter by missingness
                nb_missing = len(missing_pattern.findall(records[i].id))
                nb_samples = len(sample_pattern.findall(records[i].id))
                missing_ratio = nb_missing / nb_samples
                if missing_ratio >= max_miss:
                    i += 2
                    continue
                    
                for sample in sample_pattern.finditer(records[i].id):
                    if cons == True:
                        ID, record_id, record_seq = consense_fasta(records[i], records[i+1], sample, 'c', loc)
                        output_handle.write(''.join(record_id + '\n'))
                        output_handle.write(''.join(record_seq + '\n'))
                        s[r.get(ID)] = ''.join(s.get(r.get(ID)) + record_seq)
                        continue
                        
                    if sample.group(2) == '0':
                        ID, record_id, record_seq = extract_fasta(records[i], sample, '0', loc)
                    elif sample.group(2) == '1':
                        ID, record_id, record_seq = extract_fasta(records[i+1], sample, '0', loc)
                    else:
                        ID, record_id, record_seq = extract_fasta(records[i], sample, '0', loc)
                        
                    # write fasta to a common file
                    output_handle.write(''.join(record_id + '\n'))
                    output_handle.write(''.join(record_seq + '\n'))
                    s0[r.get(ID)] = ''.join(s0.get(r.get(ID)) + record_seq)
                        
                    if sample.group(3) == '0':
                        ID, record_id, record_seq = extract_fasta(records[i], sample, '1', loc)
                    elif sample.group(3) == '1':
                        ID, record_id, record_seq = extract_fasta(records[i+1], sample, '1', loc)
                    else:
                        ID, record_id, record_seq = extract_fasta(records[i+1], sample, '1', loc)
                        
                    # write fasta to a common file
                    output_handle.write(''.join(record_id + '\n'))
                    output_handle.write(''.join(record_seq + '\n'))
                    s1[r.get(ID)] = ''.join(s1.get(r.get(ID)) + record_seq)
                
                output_handle.write('\n')
                i += 2
                j += 1
                partition = 'CHARSET ' + 'p' +  str(j) + '=' + str(l+1) + '-' + str(l+len(record_seq)) + ';'
                partition_handle.write(''.join(partition + '\n'))
                l += len(record_seq)
            
# write out fasta format
if cons == True:
    make_fasta(s, path, fasta_out)
else:
    make_fasta_alleles(s0, s1, path, fasta_out)

# write out nexus format
if cons == True:
    make_nexus(s, path, fasta_out)
else:
    make_nexus_alleles(s0, s1, path, fasta_out)

# write out phylip format
if cons == True:
    make_phylip(s, path, fasta_out)
else:
    make_phylip_alleles(s0, s1, path, fasta_out)

print("********".format())
print("Number of loci available: {0}".format(int(i/2)))
print("Number of loci extracted: {0}".format(j))

