import sys
import numpy as np
from blosum62 import blosum62

# arguments
SEQ_TYPE = "dna"
SEQ1 = ""
SEQ2 = ""
MATCH_VAL = 1
MISMATCH_VAL = -1
GAP_VAL = -2


def main():
    # check if arguments are empty
    if len(sys.argv) == 1:
        print_help("There were no arguments passed.")
        exit(1)
    # check for help screen
    if sys.argv[1] == "--help":
        print_help()
        exit(0)
    # make sure len arguments is even (filename counts as 1st so test if even)
    if len(sys.argv) % 2 == 0:
        print_help(
            "Your argument length is not valid. Please re-check your passed arguments"
        )
        exit(1)
    # read arguments
    for i in range(1, len(sys.argv), 2):
        handle_arg(sys.argv[i], sys.argv[i + 1])

    # run global sequence alignment

    # create empty score and traceback array
    scores = np.zeros((len(SEQ2) + 1, len(SEQ1) + 1))
    traceback = np.zeros((len(SEQ2) + 1, len(SEQ1) + 1), dtype="S4")

    # init first column and row with default values
    for i in range(len(scores[0])):
        scores[0][i] = i * GAP_VAL
        traceback[0][i] = "left"
    for i in range(len(scores.T[0])):
        scores.T[0][i] = i * GAP_VAL
        traceback.T[0][i] = "up"

    traceback[0][0] = "done"

    # run forward pass
    # loop through the remaining columns and rows and calculate the values
    for row in range(1, len(scores)):
        for col in range(1, len(scores.T)):
            scores[row][col], traceback[row][col] = handle_idx(col, row, scores)

    # traceback
    seq1 = ""
    seq2 = ""
    col = len(scores.T) - 1
    row = len(scores) - 1
    current_dir = ""

    while True:
        # get the direction of the current cell
        current_dir = traceback[row][col].decode("UTF-8")
        if current_dir == "diag":
            # grab each seq value and increment both rows
            seq1 = SEQ1[col - 1] + seq1
            seq2 = SEQ2[row - 1] + seq2
            row -= 1
            col -= 1
        elif current_dir == "up":
            # gap for seq1
            seq1 = "-" + seq1
            seq2 = SEQ2[row - 1] + seq2
            row -= 1
        elif current_dir == "left":
            # gap for seq2
            seq1 = SEQ1[col - 1] + seq1
            seq2 = "-" + seq2
            col -= 1
        elif current_dir == "done":
            # back at index
            break
        else:
            # should not get here
            print("ERROR: There was an unknown error")
            exit(1)

    # calulate the score
    score = 0
    for i in range(len(seq1)):
        # impose penalty for a gap
        if seq1[i] == "-" or seq2[i] == "-":
            score += GAP_VAL
        else:
            if SEQ_TYPE == "dna":
                if seq1[i] == seq2[i]:
                    # match
                    score += MATCH_VAL
                else:
                    # mismatch
                    score += MISMATCH_VAL
            else:
                # get blosum value from map
                score += blosum62[f"{seq1[i]}{seq2[i]}"]

    print("--Global Sequence Alignment in Python--")
    print(f"Sequencing Type: {SEQ_TYPE}")
    print(f"Sequence 1: {SEQ1}")
    print(f"Sequence 2: {SEQ2}")
    print(f"Match: {MATCH_VAL}")
    print(f"Mismatch: {MISMATCH_VAL}")
    print(f"Gap Penalty: {GAP_VAL}")
    print("--OUTPUT--")
    print(seq1)
    print(seq2)
    print(f"Total Score: {score}")


def handle_idx(col: int, row: int, vec):
    # calculate value for left, top, and diagonal
    left = vec[row][col - 1] + GAP_VAL
    top = vec[row - 1][col] + GAP_VAL
    diag = vec[row - 1][col - 1]

    # check for matches
    if SEQ_TYPE == "dna":
        # check for match in diag
        if SEQ1[col - 1] == SEQ2[row - 1]:
            diag += MATCH_VAL
        else:
            diag += MISMATCH_VAL
    else:
        # get blosum62 scoring value
        diag += blosum62[f"{SEQ1[col - 1]}{SEQ2[row - 1]}"]

    # check for max value
    if left > top and left > diag:
        return left, "left"
    elif top > left and top > diag:
        return top, "up"
    elif diag > top and diag > left:
        return diag, "diag"
    else:
        # there are equal values, determine best one to return
        if diag >= top and diag >= left:
            return diag, "diag"
        elif left > top and left > diag:
            return left, "left"
        else:
            return top, "up"


def handle_arg(arg1: str, arg2: str):
    # inherit global scope
    global SEQ_TYPE
    global SEQ1
    global SEQ2
    global MATCH_VAL
    global MISMATCH_VAL
    global GAP_VAL

    # clean inouts
    arg1 = arg1.strip()
    arg2 = arg2.strip()

    # parse the command line args
    if arg1 == "--seqtype":
        if arg2 != "dna" and arg2 != "protein":
            print_help(
                f"--seqtype only accepts the values (dna) or (protein). Passed: ({arg2})"
            )
            exit(1)
        SEQ_TYPE = arg2
        if SEQ_TYPE == "protein":
            GAP_VAL = -10
    elif arg1 == "--seq1":
        SEQ1 = arg2
    elif arg1 == "--seq2":
        SEQ2 = arg2
    elif arg1 == "--match":
        MATCH_VAL = int(arg2)
    elif arg1 == "--mismatch":
        MISMATCH_VAL = int(arg2)
    elif arg1 == "--gap":
        GAP_VAL = int(arg2)


def print_help(err: str = None):
    if err:
        print("ERROR: " + err + "\n")
    print("--Global Sequence Alignment in Python--")
    print("Available Arguments:")
    print("--seqtype (str): `dna` or `protein` |DEFAULT: dna|")
    print("--seq1 (str): First sequence |REQUIRED|")
    print("--seq2 (str): Second sequence |REQUIRED|")
    print("--match (int) (dna only): Score for a match |DEFAULT: 1|")
    print("--mismatch (int) (dna only): Penalty for a mismatch |DEFAULT: -1|")
    print("--gap (int): Penalty for a gap |DEFAULT: (dna: -2) (protein: -10)|")
    print("--help: Print this help screen")


if __name__ == "__main__":
    main()
