import pandas as pd
from Bio import SeqIO
import re
import numpy as np
import seaborn as sns
from fuzzysearch import find_near_matches
import sys

def extract_sequence(dna_sequence, start_pattern, end_pattern,mismatch=1):
    # Find start pattern with one mutation
    start_matches = find_near_matches(start_pattern, dna_sequence, max_l_dist=mismatch)
    if not start_matches:
        return None
    
    # Find end pattern with one mutation
    end_matches = find_near_matches(end_pattern, dna_sequence, max_l_dist=mismatch)
    if not end_matches:
        return None
    
    # Find the closest start and end matches
    start_match = min(start_matches, key=lambda match: match.start)
    end_match = min(end_matches, key=lambda match: match.start)
    
    # Extract the sequence between start and end matches
    extracted_sequence = dna_sequence[start_match.start:end_match.start]
    
    return extracted_sequence,start_match.start

def fuzzy_match(df,prefix,start_pattern,end_pattern):
    fuzzy_seq=[]
    fuzzy_start_index=[]
    for seq in df.seq:
        result=extract_sequence(seq,start_pattern,end_pattern)
        if result is not None:
            fuzzy_seq.append(result[0])
            fuzzy_start_index.append(result[1])
        else:
            fuzzy_seq.append('')
            fuzzy_start_index.append(-1)
    df[f'{prefix}_fuzzy_seq']=fuzzy_seq
    df[f'{prefix}_fuzzy_start_index']=fuzzy_start_index
    return df

def convert_minus_to_plus(dna_sequence):
    complementary_bases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    plus_sequence = ''.join(complementary_bases.get(base, base) for base in dna_sequence)
    return plus_sequence[::-1]

def translate_dna_to_aa(dna_sequence):
    statandard_genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    genetic_code={
        'GCT':'A','GCC':'A','GCA':'A','GCG':'A', 'GCN':'A',
        'TGT':'C','TGC':'C', 'TGY':'C',
        'GAT':'D','GAC':'D','GAY':'D',
        'GAA':'E','GAG':'E','GAR':'E',
        'TTT':'F','TTC':'F', 'TTY':'F',
        'GGT':'G','GGC':'G','GGA':'G','GGG':'G', 'GGN':'G',
        'CAT':'H','CAC':'H', 'CAY':'H',
        'ATT':'I','ATC':'I','ATA':'I', 'ATH':'I',
        'AAA':'K','AAG':'K', 'AAR':'K',
        'TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L', 'YTR':'L', 'CTN':'L',
        'ATG':'M',
        'AAT':'N','AAC':'N', 'AAY':'N',
        'CCT':'P','CCC':'P','CCA':'P','CCG':'P', 'CCN':'P',
        'CAA':'Q','CAG':'Q', 'CAR':'Q',
        'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R', 'MGR':'R', 'CGN':'R', 
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S', 'TCN':'S', 'AGY':'S',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T', 'ACN':'T',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V', 'GTN':'V',
        'TGG':'W',
        'TAT':'Y','TAC':'Y', 'TAY':'Y',
        'TAA':'*','TAG':'*','TGA':'*', 'TRA':'*',
        
    }
    if dna_sequence==np.nan:
        return np.nan
    if dna_sequence=="":
        return ""
    if len(dna_sequence) % 3 != 0:
        return "length problem"
    if dna_sequence[0:3]!='ATG':
        return "start problem"
    codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]
    amino_acids = []
    for codon in codons:
        amino_acid = genetic_code.get(codon, 'X')
        # if amino_acid == '*':
        #     break  # Stop translation at the first stop codon
        amino_acids.append(amino_acid)
    protein_sequence = ''.join(amino_acids)
    return protein_sequence

## calculate the different sites between sequences
def mut_site_num(seq1,seq2):
    if seq1=="" or seq2=="":
        return np.nan
    if len(seq1)!=len(seq2):
        return np.nan
    return sum([1 for i in range(len(seq1)) if seq1[i]!=seq2[i]])

def decide_fwd_reverse(seq,fwd_start,rev_start):
    if seq[0]==fwd_start[0]:
        return 'fwd'
    elif seq[0]==rev_start[0]:
        return 'rev'
    elif seq[1]==fwd_start[1]:
        return 'fwd'
    elif seq[1]==rev_start[1]:
        return 'rev'
    else:
        return "-1"