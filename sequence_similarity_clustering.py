from Bio import pairwise2
import numpy as np
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage, single
import random
import pandas as pd

def get_normalized_smith_waterman_score(s1,s2):
    score12=pairwise2.align.localxx(s1,s2,score_only=1)
    score1=pairwise2.align.localxx(s1,s1,score_only=1)
    score2=pairwise2.align.localxx(s2,s2,score_only=1)
    return score12/np.sqrt(score1*score2)

def get_cluster_idx(matrix,threshold):
    P_dist = []
    for i in range(matrix.shape[0]):
        P_dist += (1-matrix[i,(i+1):]).tolist()
    P_dist=np.array(P_dist)
    P_link=single(P_dist)
    cluster=fcluster(P_link,threshold,'distance')
    print('Number of clusters:',max(cluster))
    return cluster

def cluster_idx_seq_list(unique_seq1,threshold):
    matrix=np.array([[get_normalized_smith_waterman_score(s1,s2) for s2 in unique_seq1] for s1 in unique_seq1])
    cluster_idx=get_cluster_idx(matrix,threshold)
    return cluster_idx


df=pd.read_csv('data/seq.csv')
unique_seq1=df.seq1.unique().tolist()
cluster_idx1=cluster_idx_seq_list(unique_seq1,threshold=0.3)