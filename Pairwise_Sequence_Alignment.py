import sys
from GlobalAlignment import global_alignment
from Semi_Global_Alignment import semi_global_alignment
from Local_Alignment import local_alignment

def inputValidation(sortedParameters):
    if len(sys.argv) != 11:
        return False
    else:
        "The checks for flags in input"
        for i in range(1,10,2):
            if not(flag(sys.argv[i],i,sortedParameters)):
                print(sys.argv[i])
                return False
        for j in range(0,5,1):
            if sortedParameters[j] is None:
                return False
    return True


def flag(argument,index, sP):
    if argument == "-p":
        if sys.argv[index+1] == "T":
            sP[2] = "T"
            return True
        elif sys.argv[index+1] == "F":
            sP[2] = "F"
            return True
        else:
            return False
    elif argument == "-i":
        sP[0] = sys.argv[index+1]
        return True
    elif argument == "-j":
        Seq2 = True
        sP[1] = sys.argv[index + 1]
        return True
    elif argument == "-atype":
        if sys.argv[index + 1] == "G":
            sP[3] = sys.argv[index + 1]
            return True
        elif sys.argv[index + 1] == "S":
            sP[3] = sys.argv[index + 1]
            return True
        elif sys.argv[index + 1] == "L":
            sP[3] = sys.argv[index + 1]
            return True
        else:
            print("Incorrect Atype")
            return False

    elif argument == "-o":
        sP[4] = sys.argv[index + 1]
        return True
    else:
        print("Invalid flag or argument passed valid flag")
        return False


def sequenceToString(sequence_File_Name,SequenceType):
    sequence_String = ""
    if SequenceType == "T":
        try:
            sequenceFile = open("Protein_Sequences/"+sequence_File_Name, "r")
            """removing the header info"""
            sequence_String = sequenceFile.readline()
            sequence_String = ""
            for line in sequenceFile:
                sequence_String += line.strip("\n")

        except IOError:
            print("Could not read file:", sequence_File_Name)
            sys.exit()

    elif SequenceType == "F":
        try:

            sequenceFile = open("DNA_Sequences/"+sequence_File_Name, "r")
            "removing the header info"
            sequence_String = sequenceFile.readline()
            sequence_String = ""
            for line in sequenceFile:
                sequence_String += line.strip("\n")

        except IOError:
            print("Could not read file:", sequence_File_Name)
            sys.exit()
    else:
        print("Could not read file:", sequence_File_Name)
        sys.exit()

    return sequence_String


def getScoringMatrix(matrixType, scoringMatrix):
    if matrixType == "T":
        try:
            matrix_file = open("Scoring_Matricies/BLOSUM45", "r")
            for i in range(0,25,1):
                scoringMatrix.append([0]*25)
            gapScore = getMatrix(scoringMatrix, matrix_file)
        except IOError:
            print("Could not open Scoring Matrix")
            sys.exit()
    elif matrixType == "F":
        try:
            matrix_file = open("Scoring_Matricies/dnaMatrix", "r")
            for i in range(0,5,1):
                scoringMatrix.append([0]*5)
            gapScore = getMatrix(scoringMatrix,matrix_file)
        except IOError:
            print("Could not open Scoring Matrix")
            sys.exit()
    else:
        print("Check parameters and try again")
        sys.exit()
    return gapScore


def getMatrix(scoringMatrix,matrix_file):
    firstLine = True
    row = 0
    for line in matrix_file:
        # print("currentLine: ",line.split(" "))

        if line[0] != "#":
            if len(line) < 5:
                gapScore = line.strip("\n")
            elif firstLine == True:
                count = 1
                for i in range(0, len(line), 1):
                    if line[i].isalnum():
                        scoringMatrix[row][count] = line[i]
                        count += 1
                firstLine = False
                #print(scoringMatrix)

            else:
                currRow = line.strip("\n")
                currRow = currRow.strip('').split(" ")
                currRow = [x for x in currRow if x != ""]
                scoringMatrix[row] = currRow
            row += 1

    return gapScore



def generate_alignmentMatrix(long_seq, seq_two):

    l_seq_size_padded = len(long_seq) + 2 # the additional 2 is required for 0 and gap char
    s_seq_size_padded = len(seq_two) + 2
    alingment_matrix = [0]*s_seq_size_padded

    for col in range(0,s_seq_size_padded):
        alingment_matrix[col] = [0]*l_seq_size_padded
    for row in range(0,len(seq_two),1):
        alingment_matrix[row+2][0] = seq_two[row]

    alingment_matrix[0] = [0,"-"] + list(long_seq)
    alingment_matrix[1][0] = "-"
    return alingment_matrix

def perform_alignment(Sequence_0, Sequence_1, ScoringMatrix,alignment_type,gap_score ):
    if alignment_type == "G" or alignment_type == "S" or alignment_type == "L":
        longest_sequence = Sequence_0
        sequnce_two = Sequence_1
        if len(Sequence_1) > len(Sequence_0):
            longest_sequence = Sequence_1
            sequnce_two = Sequence_0

        alingment_matrix = generate_alignmentMatrix(longest_sequence,sequnce_two)

        if alignment_type == "G":
            print("Global Alignment")
            global_alignment(longest_sequence, sequnce_two,alingment_matrix,ScoringMatrix,gap_score)
        """
        elif alignment_type == "S":
            print("Semi-Global Alignment")
            semi_global_alignment(longest_sequence, sequnce_two,alingment_matrix,ScoringMatrix)
        elif alignment_type == "L":
            print("Local Alignemnt")
            local_alignment(longest_sequence, sequnce_two,alingment_matrix,ScoringMatrix)
        else:
            print("Alignment type parameter is incorrect")
            sys.exit()
            """
    return "TACO"

def printScoringMatrix(scoringMatrix):
    for row in range(0,len(scoringMatrix),1):
        print(scoringMatrix[row],"\tlen: ",len(scoringMatrix[row]))


def main():
    sortedParameters = [None] * 5
    if inputValidation(sortedParameters):
        Sequence_0_FName = sortedParameters[0]
        Sequence_1_FName = sortedParameters[1]
        Sequence_0 = sequenceToString(Sequence_0_FName, sortedParameters[2])
        Sequence_1 = sequenceToString(Sequence_1_FName, sortedParameters[2])
        ScoringMatrix = []
        gap_score = getScoringMatrix(sortedParameters[2],ScoringMatrix)
        info = perform_alignment(Sequence_0, Sequence_1, ScoringMatrix,sortedParameters[3],gap_score )

    else:
        print("Parameter Error")
        sys.exit()

main()