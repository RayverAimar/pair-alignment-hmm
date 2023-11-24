import numpy as np

def pair_alignment_hmm(sequence1, sequence2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    transition_probs = {
        'M': {'M': 0.8, 'X': 0.1, 'Y': 0.1},
        'X': {'M': 0.4, 'X': 0.4, 'Y': 0.2},
        'Y': {'M': 0.4, 'X': 0.2, 'Y': 0.4}
    }
    
    initial_probs = {'M': 1.0, 'X': 0.0, 'Y': 0.0}
    
    score_matrix = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
    
    viterbi_matrix = np.zeros((len(sequence1) + 1, len(sequence2) + 1), dtype=int)
    
    state = 'M'
    
    for i in range(1, len(sequence1) + 1):
        for j in range(1, len(sequence2) + 1):
            scores = [score_matrix[i-1, j-1] + (match_score if sequence1[i-1] == sequence2[j-1] else mismatch_score),
                      score_matrix[i-1, j] + gap_penalty,
                      score_matrix[i, j-1] + gap_penalty]
            
            score_matrix[i, j] = max(scores)
            viterbi_matrix[i, j] = np.argmax(scores)
    
    alignment1 = ''
    alignment2 = ''
    
    i, j = len(sequence1), len(sequence2)
    
    while i > 0 or j > 0:
        if viterbi_matrix[i, j] == 0:  # Match or mismatch
            alignment1 = sequence1[i-1] + alignment1
            alignment2 = sequence2[j-1] + alignment2
            i -= 1
            j -= 1
        elif viterbi_matrix[i, j] == 1:  # Gap in sequence1
            alignment1 = sequence1[i-1] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        else:  # Gap in sequence2
            alignment1 = '-' + alignment1
            alignment2 = sequence2[j-1] + alignment2
            j -= 1
    
    return alignment1, alignment2

sequence1 = "MAFTWGLLALM"
sequence2 = "MATWGLAM"
alignment1, alignment2 = pair_alignment_hmm(sequence1, sequence2)
print("Secuencia 1:", alignment1)
print("Secuencia 2:", alignment2)