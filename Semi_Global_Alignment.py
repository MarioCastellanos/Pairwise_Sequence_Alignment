def printScoringMatrix(scoringMatrix):
    for row in range(0, len(scoringMatrix), 1):
        print(scoringMatrix[row], "\tlen: ", len(scoringMatrix[row]))


def semi_global_alignment(alignment_matrix,ScoringMatrix,trace_matrix,gap_score):
    print("asdf")
    printScoringMatrix(alignment_matrix)

    #fill_matrix_score(alignment_matrix, trace_matrix, scoring_matrix, gap_score)
