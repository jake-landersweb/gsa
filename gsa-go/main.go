package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"
)

// arguments
var SEQ_TYPE string = "dna"
var SEQ1 string = ""
var SEQ2 string = ""
var MATCH_VAL float64 = 1
var MISMATCH_VAL float64 = -1
var GAP_VAL float64 = -2
var GAP_EXTENSION float64 = -0.5
var OUTPUT_TYPE string = "verbose"

func main() {
	// check if arguments are empty
	if len(os.Args) == 1 {
		print_help("There were no arguments passed", 1)
	}
	// check for help screen
	if os.Args[1] == "--help" {
		print_help("", 0)
	}
	// make sure length is even
	if len(os.Args)%2 == 0 {
		print_help("Your argument length is not valid. Please re-check your passed arguments", 1)
	}
	// read arguments
	idx := 1
	for {
		handle_arg(os.Args[idx], os.Args[idx+1])
		idx += 2
		if idx >= len(os.Args) {
			break
		}
	}

	// create empty arrays
	scores := make([][]float64, len(SEQ2)+1)
	traceback := make([][]string, len(SEQ2)+1)
	for i := range scores {
		scores[i] = make([]float64, len(SEQ1)+1)
		traceback[i] = make([]string, len(SEQ1)+1)
	}

	// init first column and row with default values
	for i := range scores[0] {
		scores[0][i] = float64(i) * GAP_VAL
		traceback[0][i] = "left"
	}
	for i := range scores {
		scores[i][0] = float64(i) * GAP_VAL
		traceback[i][0] = "up"
	}
	traceback[0][0] = "done"

	// forward pass
	// loop through remaining rows and calculate the values
	for row := 1; row < len(scores); row++ {
		for col := 1; col < len(scores[0]); col++ {
			scores[row][col], traceback[row][col] = handle_idx(col, row, &scores, &traceback)
		}
	}

	// traceback
	seq1 := ""
	seq2 := ""
	col := len(scores[0]) - 1
	row := len(scores) - 1
	current_dir := ""

	for {
		// get the direction of the current cell
		current_dir = traceback[row][col]
		if current_dir == "diag" {
			// grab each seq value and increment both rows
			seq1 = fmt.Sprintf("%v%v", string(SEQ1[col-1]), seq1)
			seq2 = fmt.Sprintf("%v%v", string(SEQ2[row-1]), seq2)
			row -= 1
			col -= 1
		} else if current_dir == "up" {
			// gap for seq1
			seq1 = fmt.Sprintf("-%v", seq1)
			seq2 = fmt.Sprintf("%v%v", string(SEQ2[row-1]), seq2)
			row -= 1
		} else if current_dir == "left" {
			// gap for seq2
			seq1 = fmt.Sprintf("%v%v", string(SEQ1[col-1]), seq1)
			seq2 = fmt.Sprintf("-%v", seq2)
			col -= 1
		} else if current_dir == "done" {
			// back at index
			break
		} else {
			// should not get here
			fmt.Println("ERROR: There was an unknown error")
			os.Exit(1)
		}
	}

	// calulate the score
	score := 0.0
	for i := 0; i < len(seq1); i++ {
		// impose penalty for a gap
		if string(seq1[i]) == "-" || string(seq2[i]) == "-" {
			if i == 0 {
				// no possible extension
				score += GAP_VAL
			} else {
				// more restrictive way to calculate gap extension
				// if (string(seq1[i]) == "-" && string(seq1[i-1]) == "-") || (string(seq2[i]) == "-" && string(seq2[i-1]) == "-") {

				// check if there is a gap extension from either strand
				if string(seq1[i]) == "-" || string(seq2[i]) == "-" || string(seq1[i-1]) == "-" || string(seq2[i-1]) == "-" {
					score += GAP_EXTENSION
				} else {
					score += GAP_VAL
				}
			}
		} else {
			if SEQ_TYPE == "dna" {
				if seq1[i] == seq2[i] {
					// match
					score += MATCH_VAL
				} else {
					// mismatch
					score += MISMATCH_VAL
				}
			} else {
				// get blosum value from map
				score += blosum62[fmt.Sprintf("%v%v", string(seq1[i]), string(seq2[i]))]
			}
		}

	}

	// print results
	if OUTPUT_TYPE == "verbose" {
		fmt.Println("--Global Sequence Alignment in Go--")
		fmt.Println(fmt.Sprintf("Sequence Type: %s", SEQ_TYPE))
		fmt.Println(fmt.Sprintf("Sequence 1   : %s", SEQ1))
		fmt.Println(fmt.Sprintf("Sequence 2   : %s", SEQ2))
		fmt.Println(fmt.Sprintf("Match        : %v", MATCH_VAL))
		fmt.Println(fmt.Sprintf("Mismatch     : %v", MISMATCH_VAL))
		fmt.Println(fmt.Sprintf("Gap Penalty  : %v", GAP_VAL))
		fmt.Println(fmt.Sprintf("Gap Extension: %v", GAP_EXTENSION))
		fmt.Println("--OUTPUT--")
		print_alignment(&seq1, &seq2)
		fmt.Println(fmt.Sprintf("Total Score: %v", score))
	} else if OUTPUT_TYPE == "onlyscore" {
		fmt.Println(fmt.Sprintf("Total Score: %v", score))
	}
}

func print_alignment(seq1 *string, seq2 *string) {
	// grab constants for list manip
	length := len(*seq1)
	beginning := 0
	end := 50

	// num base pairs iterated over for each strand
	char1 := 1
	char2 := 1

	// print 50 lines of each seq at a time
	for {
		// create temp lists
		tmp1 := (*seq1)[beginning:min(length, end)]
		tmp2 := (*seq2)[beginning:min(length, end)]

		// compose lines
		line1 := ""
		line2 := ""
		line3 := ""

		for i := range tmp1 {
			line1 += string(tmp1[i])
			line2 += string(tmp2[i])
			// check for exact matches
			if tmp1[i] == tmp2[i] {
				line3 += "*"
			} else {
				line3 += " "
			}
			// only add char if not a gap
			if string(tmp1[i]) != "-" {
				char1 += 1
			}
			if string(tmp2[i]) != "-" {
				char2 += 1
			}
		}
		// add traversed chars after each sequence line
		line1 += fmt.Sprintf(" |%v", char1)
		line2 += fmt.Sprintf(" |%v", char2)
		fmt.Println(line1)
		fmt.Println(line2)
		fmt.Println(line3)
		// increase list indexes
		beginning += 50
		end += 50
		if beginning > length {
			break
		}
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func handle_idx(col int, row int, scores *[][]float64, traceback *[][]string) (float64, string) {
	// for calculating without extension
	// left := (*scores)[row][col-1] + GAP_VAL
	// top := (*scores)[row-1][col] + GAP_VAL
	// diag := (*scores)[row-1][col-1]

	// calculate with gap extension
	left := (*scores)[row][col-1]
	top := (*scores)[row-1][col]
	if (*traceback)[row][col] != "diag" {
		left += GAP_EXTENSION
		top += GAP_EXTENSION
	} else {
		left += GAP_VAL
		top += GAP_VAL
	}
	diag := (*scores)[row-1][col-1]

	if SEQ_TYPE == "dna" {
		// check for match in diag
		if SEQ1[col-1] == SEQ2[row-1] {
			diag += MATCH_VAL
		} else {
			diag += MISMATCH_VAL
		}
	} else {
		// add blosum62 value from map
		diag += blosum62[fmt.Sprintf("%v%v", string(SEQ1[col-1]), string(SEQ2[row-1]))]
	}

	// check for max value
	if left > top && left > diag {
		return left, "left"
	} else if top > left && top > diag {
		return top, "up"
	} else if diag > top && diag > left {
		return diag, "diag"
	} else {
		// there are equal values, determine best one to return
		if diag >= top && diag >= left {
			return diag, "diag"
		} else if left > top && left > diag {
			return left, "left"
		} else {
			return top, "up"
		}
	}
}

func handle_arg(arg1 string, arg2 string) {
	// clean the arguments
	arg1 = strings.TrimSpace(arg1)
	arg2 = strings.TrimSpace(arg2)

	// parse the command line args
	if arg1 == "--seqtype" {
		if arg2 != "dna" && arg2 != "protein" {
			print_help(fmt.Sprintf("--seqtype only accepts the values (dna) or (protein). Passed: (%s)", arg2), 1)
		}
		SEQ_TYPE = arg2
		if SEQ_TYPE == "protein" {
			GAP_VAL = -10
		} else {
			GAP_EXTENSION = GAP_VAL
		}
	} else if arg1 == "--seq1" {
		// remove last digit if it is a *
		if last := len(arg2) - 1; last >= 0 && arg2[last] == '*' {
			arg2 = arg2[:last]
		}
		SEQ1 = strings.ToUpper(arg2)
	} else if arg1 == "--seq2" {
		// remove last digit if it is a *
		if last := len(arg2) - 1; last >= 0 && arg2[last] == '*' {
			arg2 = arg2[:last]
		}
		SEQ2 = strings.ToUpper(arg2)

	} else if arg1 == "--match" {
		tmp, err := strconv.ParseFloat(arg2, 64)
		if err != nil {
			print_help(err.Error(), 1)
		}
		MATCH_VAL = tmp
	} else if arg1 == "--mismatch" {
		tmp, err := strconv.ParseFloat(arg2, 64)
		if err != nil {
			print_help(err.Error(), 1)
		}
		MISMATCH_VAL = -tmp
	} else if arg1 == "--gap" {
		tmp, err := strconv.ParseFloat(arg2, 64)
		if err != nil {
			print_help(err.Error(), 1)
		}
		GAP_VAL = -tmp
	} else if arg1 == "--gapext" {
		tmp, err := strconv.ParseFloat(arg2, 64)
		if err != nil {
			print_help(err.Error(), 1)
		}
		GAP_EXTENSION = -tmp
	} else if arg1 == "--outtype" {
		if arg2 != "verbose" && arg2 != "onlyscore" {
			print_help(fmt.Sprintf("--outtype only accepts the values (verbose) or (onlyscore). Passed: (%s)", arg2), 1)
		}
		OUTPUT_TYPE = arg2
	}
}

func print_help(err string, code int) {
	if err != "" {
		fmt.Println(fmt.Sprintf("ERROR: %s", err))
	}
	fmt.Println("--Global Sequence Alignment in Go--")
	fmt.Println("Available Arguments:")
	fmt.Println("--seqtype (str): `dna` or `protein` |DEFAULT: dna|")
	fmt.Println("--seq1 (str): First sequence |REQUIRED|")
	fmt.Println("--seq2 (str): Second sequence |REQUIRED|")
	fmt.Println("--match (float) (dna only): Score for a match |DEFAULT: 1|")
	fmt.Println("--mismatch (float) (dna only): Penalty for a mismatch |DEFAULT: -1|")
	fmt.Println("--gap (float): Penalty for a gap |DEFAULT: (dna: -2) (protein: -10)|")
	fmt.Println("--gapext (float): Penalty for a gap extension |DEFAULT: (dna: -2) (protein: -0.5)|")
	fmt.Println("--outtype (str): Style of output, `verbose` or `compact` |DEFAULT: verbose|")
	fmt.Println("--help: fmt.Println this help screen")
	os.Exit(code)
}
