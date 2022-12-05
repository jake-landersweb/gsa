# Global Sequence Alignment Algorithm

Implementation of the [Global Sequence Alignment algorithm](https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm) in various programming languages.

## Python (Human):
- Built to prototype various ideas for the implementation

Run `python main.py --help` to see available arguments and how to run this program.

Example: `python main.py --seq1 $seq1 --seq2 $seq2 --seqtype dna`

## Go:
- Same features as Python +:
- Better printing
- Gap extension penalties

Make sure you have Go installed and in your path. You can download it from [https://go.dev/doc/install](https://go.dev/doc/install).

Build: `cd gsa-go && ./build.sh`

Run `./gsa-go --help` to see available arguments and how to run this program.

Example: `./gsa-go --seq1 $seq1 --seq2 $seq2 --seqtype dna`

## Python (AI)

After finishing my implementation, I decided to use [chat.openai.com](chat.openai.com) to see if it could complete this assignment for me by passing in the instructions and asking it to complete this for me. It was able to write it hauntingly well. This can be found in `gpt.py`.

The AI was able to generate the complete function that handles both DNA and proteins, and an entire BLOSUM62 substitution matrix that is uses itself for the calculations. The only piece of code I wrote was a function to give me better printing of the sequences.

Here is the transcript from the chat ai:

Me:
> This project will entail implementing the Global Sequence Alignment dynamic programming algorithm in Python (w/o using modules containing these algorithms). The program should be able to align DNA and Protein sequences. Protein sequence alignments should be scored using BLOSUM62. The programs should be able to run using Python 3 (e.g., Python 3.10, etc.).

Bot:
> To implement the Global Sequence Alignment algorithm in Python, you will need to first understand how the algorithm works. Global Sequence Alignment is a method for aligning two sequences of nucleotides or amino acids in order to determine the most similar regions between them. This is often used in the field of bioinformatics to compare DNA or protein sequences and identify similarities or differences between them.

> To implement the algorithm in Python, you will need to use a dynamic programming approach. This involves defining a scoring matrix that will hold the scores for each pair of aligned nucleotides or amino acids. You will then need to iterate over the two sequences, comparing each pair of nucleotides or amino acids and updating the scoring matrix with the appropriate scores.

> Once the scoring matrix has been filled in, you can then use it to trace back through the sequences and determine the optimal alignment between them. This will typically involve starting at the bottom right corner of the matrix and working your way back to the top left, following the path with the highest score.

> For protein sequence alignments, you will need to use a substitution matrix such as BLOSUM62 to determine the scores for aligning each pair of amino acids. This matrix contains pre-computed scores for aligning each pair of amino acids, which you can use to determine the scores for each pair of aligned amino acids in your sequences.

> Overall, implementing the Global Sequence Alignment algorithm in Python will require a combination of dynamic programming and knowledge of substitution matrices such as BLOSUM62. With the right approach, you should be able to implement a working version of the algorithm in Python that can be used to align DNA and protein sequences.

Me:
> Show me a potential implementation for this

Bot:
> Here is a potential implementation of the Global Sequence Alignment algorithm in Python, without using any pre-existing modules that contain the algorithm:

> gpt.py

> This implementation defines the global_alignment() function, which takes in two sequences and a substitution matrix (for protein sequence alignments) and returns the two aligned sequences. The function uses a dynamic programming approach to fill in the scoring matrix and then trace back through the matrix to determine the optimal alignment.

> The substitution matrix used in this example is for protein sequence alignments using the BLOSUM62 matrix. You can modify this matrix or use a different substitution matrix as needed for your specific use case.

> The example also includes sample sequences to align and prints the resulting alignments. You can modify these sequences as needed and use the global_alignment() function to align any sequences of your choosing.

You need to open `gpt.py` and edit the sequences, then run `python gpt.py` to run the program. The program does not support changing Match, Mismatch, or Gap scores.

## Author

Jake Landers 2022