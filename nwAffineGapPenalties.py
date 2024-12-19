'''Install this to run validation against MUSCLE'''
#pip install biopython

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
from Bio import Align

class NeedlemanWunsch:
    def __init__(self, match_score=1, mismatch_score=-1, gap_penalty=-2):
        """
        Basic Needleman-Wunsch algorithm implementation.

        Parameters:
        match_score (float): Score for matching characters
        mismatch_score (float): Score for mismatching characters
        gap_penalty (float): Penalty for gaps
        """
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty

    def _initialize_matrix(self, seq1, seq2):
        """Initialize the scoring matrix with gap penalties."""
        rows, cols = len(seq1) + 1, len(seq2) + 1
        matrix = np.zeros((rows, cols))
        for i in range(rows):
            matrix[i, 0] = i * self.gap_penalty
        for j in range(cols):
            matrix[0, j] = j * self.gap_penalty
            
        return matrix

    def _score(self, char1, char2):
        """Calculate score for matching/mismatching characters."""
        return self.match_score if char1 == char2 else self.mismatch_score

    def align(self, seq1, seq2):
        """
        Perform sequence alignment using Needleman-Wunsch algorithm.

        Parameters:
        seq1 (str): First sequence
        seq2 (str): Second sequence

        Returns:
        tuple: Aligned sequences, score, and traceback path
        """
        matrix = self._initialize_matrix(seq1, seq2)#scoring matrix
        rows, cols = len(seq1) + 1, len(seq2) + 1
        
        for i in range(1, rows):
            for j in range(1, cols):
                match_score = self._score(seq1[i-1], seq2[j-1])
                matrix[i, j] = max(
                    matrix[i-1, j-1] + match_score,
                    matrix[i-1, j] + self.gap_penalty,#gap in seq2
                    matrix[i, j-1] + self.gap_penalty#gap in seq1
                )

        #traceback the table
        aligned1, aligned2 = [], []
        i, j = rows-1, cols-1
        traceback_path = [(i, j)]
        
        while i > 0 or j > 0:
            if i > 0 and j > 0:
                scores = [
                    matrix[i-1, j-1],  
                    matrix[i-1, j],    
                    matrix[i, j-1]    
                ]
                move = np.argmax([
                    scores[0] + self._score(seq1[i-1], seq2[j-1]),
                    scores[1] + self.gap_penalty,
                    scores[2] + self.gap_penalty
                ])
            elif i > 0:
                move = 1  
            else:
                move = 2  
                
            if move == 0:  
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif move == 1:  
                aligned1.append(seq1[i-1])
                aligned2.append('-')
                i -= 1
            else:  
                aligned1.append('-')
                aligned2.append(seq2[j-1])
                j -= 1
            
            traceback_path.append((i, j))
        #reversing for traceback
        aligned1 = ''.join(reversed(aligned1))
        aligned2 = ''.join(reversed(aligned2))
        
        return aligned1, aligned2, matrix[-1, -1], traceback_path

    def visualize_alignment(self, seq1, seq2):
        """Visualize the alignment matrix using seaborn heatmap."""
        matrix = self._initialize_matrix(seq1, seq2)
        rows, cols = len(seq1) + 1, len(seq2) + 1
        
        for i in range(1, rows):
            for j in range(1, cols):
                match_score = self._score(seq1[i-1], seq2[j-1])
                matrix[i, j] = max(
                    matrix[i-1, j-1] + match_score,
                    matrix[i-1, j] + self.gap_penalty,
                    matrix[i, j-1] + self.gap_penalty
                )

        plt.figure(figsize=(10, 8))
        sns.heatmap(matrix, annot=True, fmt='.1f', cmap='coolwarm')
        plt.title('Alignment Matrix Visualization')
        plt.xlabel('Sequence 2')
        plt.ylabel('Sequence 1')
        plt.show()

def test_alignment():
    #testing simple seq
    aligner = NeedlemanWunsch(
        match_score=1,
        mismatch_score=-1,
        gap_penalty=-2
    )

    seq1 = "GCATGCU"
    seq2 = "GATTACA"

    aligned1, aligned2, score = aligner.align(seq1, seq2)

    print("Sequence 1:", seq1)
    print("Sequence 2:", seq2)
    print("\nAlignment:")
    print(aligned1)
    print(aligned2)
    print("\nAlignment Score:", score)
    aligner.visualize_alignment(seq1, seq2)

if __name__ == "__main__":
    test_alignment()