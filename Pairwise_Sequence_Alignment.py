
import sys


def inputValidation(sortedParameters):

    if len(sys.argv) != 11:
        return False
    else:
        "The checks for flags and input with falgs"
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



def main():


    sortedParameters = [None] * 5
    if inputValidation(sortedParameters):
        for param in sortedParameters:
            print(param)
    else:
        return -1


main()