from Bio.Align import substitution_matrices
from Bio import Align
import numpy as np

def get_seq_sim(seq1,seq2,aligner):
    return aligner.score(seq1,seq2)/(np.sqrt(aligner.score(seq1,seq1))*np.sqrt(aligner.score(seq2,seq2)))

def get_protein_seq_similarity(seq1, seq2, matrix_name="BLOSUM62"):
    """
    Calculate the similarity between two protein sequences
    using a given substitution matrix.
    """
    matrix = substitution_matrices.load(matrix_name)
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.mode = 'global'
    return get_seq_sim(seq1,seq2,aligner)
