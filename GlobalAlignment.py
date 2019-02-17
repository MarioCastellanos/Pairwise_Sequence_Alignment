from copy import deepcopy


def printScoringMatrix(scoringMatrix):
    for row in range(0, len(scoringMatrix), 1):
        print(scoringMatrix[row], "\tlen: ", len(scoringMatrix[row]))


def initialize_alignment_matrix_values(alignment_matrix,gap_score):
    int_gap_score = int(gap_score)

    for i in range(1, len(alignment_matrix[0])-1, 1):
        alignment_matrix[1][i+1] = i* int_gap_score
    for row in range(0, len(alignment_matrix)-1, 1):
        alignment_matrix[row+1][1] = row*int_gap_score


def max_value(row,col,a,b,alignment_matrix, scoring_matrix,gap_score):
    #print("row: ", row, "col: ", col, "a: ", a, "b: ", b)
    score_row = scoring_matrix[0].index(a)
    score_col = scoring_matrix[0].index(b)

    diagonal_val = alignment_matrix[row-1][col-1] + int(scoring_matrix[score_row][score_col])
    left_val = alignment_matrix[row][col-1] + int(gap_score)
    up_val = alignment_matrix[row-1][col] + int(gap_score)

    print("DiagnolaScore: ", diagonal_val)
    print("UpVal: ", up_val)
    print("leftVal: ",left_val)

    if diagonal_val >= left_val and diagonal_val >=up_val:
        print("D")
        return diagonal_val
    elif left_val>=diagonal_val and left_val >= up_val:
        print("L")
        return left_val
    else:
        print("U")
        return up_val


def fill_matrix_score(alignemnt_matrix, trace_matrix, scoring_matrix,gap_score):

    for row in range(2, len(alignemnt_matrix), 1):
        seq_1_val = alignemnt_matrix[row][0]
        for col in range(2,len(alignemnt_matrix[0]),1):
            seq_0_val = alignemnt_matrix[0][col]
            alignemnt_matrix[row][col] = max_value(row,col,seq_0_val,seq_1_val,alignemnt_matrix,scoring_matrix,gap_score)
            printScoringMatrix(alignemnt_matrix)
            taco = input("ENTERTOCONTINUE")


def global_alignment(sequence_1, sequence_2, alignment_matrix, scoring_matrix,gap_score):
   aligned_sequences = ["", ""]
   trace_matrix = deepcopy(alignment_matrix)
   initialize_alignment_matrix_values(alignment_matrix,gap_score)
   fill_matrix_score(alignment_matrix, trace_matrix, scoring_matrix,gap_score)
   printScoringMatrix(trace_matrix)
   printScoringMatrix(alignment_matrix)
