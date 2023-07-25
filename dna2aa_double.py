import pandas as pd
from Bio import SeqIO
import re
import numpy as np
import seaborn as sns
from fuzzysearch import find_near_matches
import sys
from dna2aa_utils import extract_sequence,fuzzy_match,convert_minus_to_plus,translate_dna_to_aa,decide_fwd_reverse
import argparse 
import os

def if_problem_aa(seq,problem_aa,start_c,end_c):
    if seq in problem_aa:
        return True
    elif end_c!=0 and seq[end_c]!="C":
        return True
    elif "*" in seq:
        return True
    elif len(seq) <= start_c+1:
        return True
    elif seq[start_c]!="C":
        return True
    else:
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq_fwd', type=str, default='../fastq_r7/fastq_gz/4_S4_L001_R1_001.fastq', help='fastq file')
    parser.add_argument('--fastq_rev', type=str, default='../fastq_r7/fastq_gz/4_S4_L001_R2_001.fastq', help='fastq file')
    parser.add_argument('--start_pattern', type=str, default='ATGGGCTGC', help='forward primer start')
    parser.add_argument('--end_pattern', type=str, default='GGGGGAGGCAGC', help='forward primer end')
    parser.add_argument('--output', type=str, default='../EME_Kras_new_count/4_S4_L001_R1_001_peptides_count.csv', help='output file')
    parser.add_argument('--start_c', type=int, default=0, help='start index of C')
    parser.add_argument('--end_c', type=int, default=-1, help='end index of C')
    parser.add_argument('--fwd_rev_mix',type=bool,default=False,help='if the fastq file is mixed with fwd and rev')
    parser.add_argument('--fwd_start_pattern', type=str, default='TA', help='forward primer start')
    parser.add_argument('--rev_start_pattern', type=str, default='GT', help='reverse primer start')
    args = parser.parse_args()
    # Read in the data
    print(args)
    input_fwd_fastq=args.fastq_fwd
    input_rev_fastq=args.fastq_rev
    start_pattern=args.start_pattern
    end_pattern=args.end_pattern
    output_file=args.output
    start_c=args.start_c
    end_c=args.end_c
    fwd_start_pattern=args.fwd_start_pattern
    rev_start_pattern=args.rev_start_pattern
    prefix=os.path.basename(input_fwd_fastq).split('.')[0]
    work_dir=os.path.dirname(output_file)
    if os.path.exists(f'{work_dir}/temp/{prefix}_id_DNA_aa.csv'):
        merged_df=pd.read_csv(f'{work_dir}/temp/{prefix}_id_DNA_aa.csv')
        print(f"Loaded {prefix}_id_DNA_aa.csv")
    
    else:
        forward_records=list(SeqIO.parse(input_fwd_fastq,'fastq'))
        reversed_records=list(SeqIO.parse(input_rev_fastq,'fastq'))
        forward_df=pd.DataFrame({'id':[record.id for record in forward_records],'seq':[str(record.seq) for record in forward_records]})
        reversed_df=pd.DataFrame({'id':[record.id for record in reversed_records],'seq':[str(record.seq) for record in reversed_records]})
        if args.fwd_rev_mix:
            print("Matching fwd and rev...")
            total_df=pd.concat([forward_df,reversed_df])
            total_df['fwd_rev']=total_df.seq.apply(lambda x: decide_fwd_reverse(x,fwd_start_pattern,rev_start_pattern))
            forward_df=total_df[total_df.fwd_rev=='fwd']
            reversed_df=total_df[total_df.fwd_rev=='rev']
            print("Didn't matched:",total_df[total_df.fwd_rev=="-1"].shape[0])

        print("Forward records: ",len(forward_records))
        print("Reversed records: ",len(reversed_records))
        reversed_df['seq']=reversed_df.seq.apply(convert_minus_to_plus)
        forward_df=fuzzy_match(forward_df,'fwd',start_pattern,end_pattern)
        reversed_df=fuzzy_match(reversed_df,'rev',start_pattern,end_pattern)
        forward_df.rename(columns={'seq':'fwd_seq'},inplace=True)
        reversed_df.rename(columns={'seq':'rev_seq'},inplace=True)
        forward_df['fwd_fuzzy_aa']=forward_df.fwd_fuzzy_seq.apply(translate_dna_to_aa)
        reversed_df['rev_fuzzy_aa']=reversed_df.rev_fuzzy_seq.apply(translate_dna_to_aa)
        print("Complete matching.")
        forward_df['fwd_matched_len']=forward_df.fwd_fuzzy_seq.apply(len)
        reversed_df['rev_matched_len']=reversed_df.rev_fuzzy_seq.apply(len)
        merged_df=pd.merge(forward_df,reversed_df, on=['id'],how='outer')
        dirs = 'temp'
        if not os.path.exists(f'{work_dir}/{dirs}'):
            os.makedirs(f'{work_dir}/{dirs}')
        merged_df.to_csv(f'{work_dir}/temp/{prefix}_id_DNA_aa.csv',index=False)
        print(f"Saved {prefix}_id_DNA_aa.csv")
    id_aa_1=[]
    id_aa_2=[]
    problem_id=[]
    problem_aa=['','start problem','length problem',np.NAN]
    for i,row in merged_df.iterrows():
        fwd_problem=if_problem_aa(row['fwd_fuzzy_aa'],problem_aa,start_c,end_c)
        rev_problem=if_problem_aa(row['rev_fuzzy_aa'],problem_aa,start_c,end_c)
        if fwd_problem and rev_problem:
            problem_id.append(row['id'])
        elif row['fwd_fuzzy_aa']==row['rev_fuzzy_aa']:
            id_aa_1.append([row['id'],row['fwd_fuzzy_aa']])
        elif fwd_problem==False and rev_problem:
            id_aa_1.append([row['id'],row['fwd_fuzzy_aa']])
        elif fwd_problem and rev_problem==False:
            id_aa_1.append([row['id'],row['rev_fuzzy_aa']])
        else:
            id_aa_2.append([row['id'],row['fwd_fuzzy_aa']])
            id_aa_2.append([row['id'],row['rev_fuzzy_aa']])
    print("Problem id: ",len(problem_id))
    print("id_aa_1: ",len(id_aa_1))
    print("id_aa_2: ",len(id_aa_2)/2)
    id_aa_1_df=pd.DataFrame(id_aa_1,columns=['id','aa'])
    id_aa_2_df=pd.DataFrame(id_aa_2,columns=['id','aa'])
    aa_count_1v1=id_aa_1_df.groupby('aa').count().reset_index()
    aa_count_1v2=id_aa_2_df.groupby('aa').count().reset_index()
    aa_count_1v1.rename(columns={'id':'count_1v1'},inplace=True)
    aa_count_1v2.rename(columns={'id':'count_1v2'},inplace=True)
    count_df=pd.merge(aa_count_1v1,aa_count_1v2,on='aa',how='outer')
    count_df['either_count']=count_df.loc[:,['count_1v1','count_1v2']].sum(axis=1,min_count=1)
    count_df['len']=count_df.aa.apply(len)
    count_df.rename(columns={'aa':'peptides'},inplace=True)
    count_df.reset_index(drop=True,inplace=True)
    count_df.to_csv(output_file,index=False)