import sys


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
                # print(scoringMatrix)

            else:
                currRow = line.strip("\n")
                currRow = currRow.strip('').split(" ")
                currRow = [x for x in currRow if x != ""]
                scoringMatrix[row] = currRow
            row += 1

    return gapScore


def printScoringMatrix(scoringMatrix):
    for row in range(0,len(scoringMatrix),1):
        print(scoringMatrix[row])


def main():
    sortedParameters = [None] * 5
    if inputValidation(sortedParameters):
        Sequence_0_FName = sortedParameters[0]
        Sequence_1_FName = sortedParameters[1]
        Sequence_0 = sequenceToString(Sequence_0_FName, sortedParameters[2])
        Sequence_1 = sequenceToString(Sequence_1_FName, sortedParameters[2])
        ScoringMatrix = []
        gap_score = getScoringMatrix(sortedParameters[2],ScoringMatrix)

    else:
        print("Parameter Error")
        sys.exit()

main()