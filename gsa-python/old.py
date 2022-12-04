def main2():
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

    # print("Running Global Sequence Alignment with the following parameters:")
    # print(f"seq_type: {SEQ_TYPE}")
    # print(f"seq1: {SEQ1}")
    # print(f"seq2: {SEQ2}")
    # print(f"match_val: {MATCH_VAL}")
    # print(f"mismatch_val: {MISMATCH_VAL}")
    # print(f"gap_val: {GAP_VAL}")

    # run global sequence alignment

    # create empty vector
    vec = np.zeros((len(SEQ2) + 1, len(SEQ1) + 1))

    # init first column and row with -2,-4,-6,...etc.
    for i in range(len(vec[0])):
        vec[0][i] = i * -2
    for i in range(len(vec.T[0])):
        vec.T[0][i] = i * -2

    # run forward pass
    # loop through the remaining columns and rows and calculate the values
    for col in range(1, len(vec.T)):
        for row in range(1, len(vec)):
            vec[row][col] = handle_idx_forward(col, row, vec)

    # run backwards to calculate sequences
    # everything is run in a list to allow for multiple paths to be created

    # init at bottom right value (score)
    cols = [len(vec.T) - 1]
    rows = [len(vec) - 1]

    # keep history of what cols and rows visited
    col_histories = [[]]
    row_histories = [[]]

    # loop backwards
    while True:
        # loop through cols and rows lists
        for i in range(len(col_histories)):
            if len(col_histories[i]) <= max(len(SEQ1), len(SEQ2)):
                # save current row and col in history
                col_histories[i].append(cols[i])
                row_histories[i].append(rows[i])

                # get the movement lists
                c, r = handle_idx_backward(cols[i], rows[i], vec)
                # check return list lengths and create new col row pair if needed
                if len(c) > 1:
                    # banching paths

                    # duplicate paths and col/row indexes n number of times
                    for _ in range(len(c) - 1):
                        cols.append(cols[i])
                        rows.append(rows[i])
                        col_histories.append(col_histories[i].copy())
                        row_histories.append(row_histories[i].copy())

                    # set the column and row index for each branching path
                    for idx in range(1, len(c) + 1):
                        cols[-idx] = c[len(c) - idx]
                        rows[-idx] = r[len(c) - idx]
                else:
                    # no branching paths
                    cols[i] = c[0]
                    rows[i] = r[0]

        # check if any more paths need to run
        if not any(
            len(col_histories[i]) <= max(len(SEQ1), len(SEQ2))
            for i in range(len(col_histories))
        ):
            break

    # compose sequences
    al1 = ["" for _ in range(len(col_histories))]
    al2 = ["" for _ in range(len(row_histories))]

    # loop through the col values
    for idx, col_history in enumerate(col_histories):
        for i in range(1, len(col_history)):
            # check for gap
            if col_history[i] == col_history[i - 1]:
                al1[idx] = "-" + al1[idx]
            else:
                al1[idx] = SEQ1[col_history[i]] + al1[idx]

    # loop through the row values
    for idx, row_history in enumerate(row_histories):
        for i in range(1, len(row_history)):
            # check for gap
            if row_history[i] == row_history[i - 1]:
                al2[idx] = "-" + al2[idx]
            else:
                al2[idx] = SEQ2[row_history[i]] + al2[idx]

    # keep track of best performing alignment
    best_a1 = ""
    best_a2 = ""
    best_score = -10000

    # calculate scores of all posible path
    for a1, a2 in zip(al1, al2):
        score = calculate_score(a1, a2)
        if score > best_score:
            best_a1 = a1
            best_a2 = a2
            best_score = score

    # print best alignment map and score
    # print(f"{best_a1} {best_a2} {best_score}")
    print("BEST:")
    print(best_a1)
    print(best_a2)
    print(best_score)


def handle_idx_forward(col: int, row: int, vec):
    # calculate value for left, top, and diagonal
    left = vec[row][col - 1] + GAP_VAL
    top = vec[row - 1][col] + GAP_VAL
    diag = vec[row - 1][col - 1]

    # check for match in diag
    if SEQ1[col - 1] == SEQ2[row - 1]:
        diag += MATCH_VAL
    else:
        diag += MISMATCH_VAL

    # return the max value
    return max([left, top, diag])


def handle_idx_backward(col: int, row: int, vec):
    # check for match
    if SEQ1[col - 1] == SEQ2[row - 1]:
        col -= 1
        row -= 1
        return [col], [row]

    # go towards highest neighbor
    top = vec[row - 1][col]
    diag = vec[row - 1][col - 1]
    left = vec[row][col - 1]

    # run logic based on whatever value is largest
    if top > diag and top > left:
        # gap
        row -= 1
        return [col], [row]
    elif diag > top and diag > left:
        # mismatch
        row -= 1
        col -= 1
        return [col], [row]
    elif left > top and left > diag:
        # gap
        col -= 1
        return [col], [row]
    else:
        # create branching lists
        if top == diag == left:
            col_1 = col
            row_1 = row - 1
            col_2 = col - 1
            row_2 = row - 1
            col_3 = col - 1
            row_3 = row
            return [col_1, col_2, col_3], [row_1, row_2, row_3]
        elif top == diag:
            col_1 = col
            row_1 = row - 1
            col_2 = col - 1
            row_2 = row - 1
            return [col_1, col_2], [row_1, row_2]
        elif top == left:
            col_1 = col
            row_1 = row - 1
            col_2 = col - 1
            row_2 = row
            return [col_1, col_2], [row_1, row_2]
        elif diag == left:
            col_1 = col - 1
            row_1 = row - 1
            col_2 = col - 1
            row_2 = row
            return [col_1, col_2], [row_1, row_2]
        else:
            print("THERE WAS A FATAL ERROR WITH YOUR ALIGNMENT")
            exit(1)


def calculate_score(al1, al2):
    score = 0
    for i in range(len(al1)):
        if al1[i] == "-" or al2[i] == "-":
            score += GAP_VAL
        else:
            if SEQ_TYPE == "dna":
                if al1[i] == al2[i]:
                    score += MATCH_VAL
                else:
                    score += MISMATCH_VAL
            else:
                score += blosum62[f"{al1[i]}{al2[i]}"]
    return score
