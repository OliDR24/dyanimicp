#!/usr/bin/python

"""Template for dynamic programming assignment.

The code in the template is compatible with both Python 2 and Python 3
When you finish this code, it should be at least compatible with Python 3.
"""

# Packages for commandline options:
import argparse
import sys
import pickle as pk

# Built-in exchange matrices.
with open('exchange_matrices/identity.pkl', 'rb') as f:
    identity = pk.load(f)

with open('exchange_matrices/pam250.pkl', 'rb') as f:
    pam250 = pk.load(f)

with open('exchange_matrices/blosum62.pkl', 'rb') as f:
    blosum62 = pk.load(f)


def get_args():
    """Collect the inputs."""
    parser = argparse.ArgumentParser(
        prog='PROG',
        usage='%(prog)s [options]',
        description='Aligning two sequences',
        epilog='The code was co-opted from Anton Feenstra\'s and'
        'modified by Cico Zhang'
    )
    parser.add_argument('-f', '--fasta', dest='fasta', metavar='FILE',
                        required=True, help='input alignment file (fasta)')
    parser.add_argument('-e,', '--exchange_matrix', dest='exchange_matrix',
                        metavar='EXCHANGE MATRIX NAME', help='Exchange '
                        'matrix: pam250, blosum62 or identity',
                        default='pam250')
    parser.add_argument('-l', '--local', dest='align_local',
                        action='store_true', help='Local alignment',
                        default=False)
    parser.add_argument('-g', '--global', dest='align_global',
                        action='store_true', help='Global alignment',
                        default=False)
    parser.add_argument('-s', '--semi_global', dest='align_semiglobal',
                        action='store_true', help='Semi-global alignment',
                        default=False)
    parser.add_argument('-p', '--penalty', dest='gap_penalty', type=int,
                        help='Gap penalty', default=2)
    parser.add_argument('-o', '--output', dest='alignment', required=True,
                        metavar='FILE', default='output.align',
                        help='The file to store the alignment')
    parser.add_argument('-m', '--score_matrix', dest='score_matrix',
                        required=True, metavar='FILE', default='output.align',
                        help='The file to store the score matrix')
    parser.add_argument('-v', dest='print_on_screen', action='store_true',
                        help='Print the output (alignment(s) and score '
                        'matrix) on the screen', default=False)

    args = parser.parse_args()

    if args.fasta is None:
        sys.exit('Error: no input file (fasta)')

    if not (args.align_local or args.align_global or args.align_semiglobal):
        sys.exit('Error: No alignment strategy is given: global, local or '
                 'semi-global')
    if args.align_local + args.align_global + args.align_semiglobal > 1:
        sys.exit('Error: More than one alignment strategy is given.')

    if args.exchange_matrix not in ['pam250', 'blosum62', 'identity']:
        sys.exit('Unknown exchange matrix ' + args.exchange_matrix)

    return args


class Sequence:
    """Stores a sequence object."""

    def __init__(self, Label="", Sequence=""):
        """Initialize a new Sequence object.

        Label -- identifier of sequence (text)
        Sequence -- sequence string in single-letter alphabet
        """
        self.Label = Label
        self.Sequence = Sequence

    # this makes that you can do 'print sequence' and get nice output:
    def __str__(self):
        """Return string representation of a Sequence object."""
        # newline-delimited values of all the attributes
        return ">%s\n%s" % (self.Label, self.Sequence)


def readSequences(lines):
    """Return Sequences object.

    lines -- list of lines or any object that behaves like it

    This routine parses a fasta file and returns a list of Sequence objects
    containing the sequences with label and sequence data set
    """
    seqs = []
    label = None
    seq_lines = []
    for line in lines:
        line = line.strip()      # strip off white space
        if not line:             # skip empty lines
            continue
        if line.startswith(';'):  # ignore comment lines
            continue
        # check for start of next sequence:
        if line.startswith('>'):  # label line
            # first, store the previous sequence if we had one:
            if seq_lines:
                seqs.append(Sequence(label, ''.join(seq_lines)))
                seq_lines = []
            # get the label (name) for the next sequence
            label = line[1:].strip()
        else:
            # collect all lines with sequence information for this sequence:
            seq_lines.append(line)
    # take care of the last sequence in the file
    seqs.append(Sequence(label, ''.join(seq_lines)))
    return seqs


def do_global_alignment(sequences, matrix, penalty):
    """Do pairwise global alignment using DP."""
  #Define the variables seq1, seq2 and the number of rows and columns for the scoring matrix
    seq1 = str(sequences[0].Sequence)
    seq2 = str(sequences[1].Sequence)
    rows = len(seq1) + 1
    cols = len(seq2) + 1


    def create_score_matrix(rows, cols, matrix, penalty):
        """Creates the score matrix using the number of rows and columns assigned previously,
        as this is a semi-global alignment initial gaps are free and as such the matrix initializes
        with zeroes in the first column and row. Again the calc_score function is used to calculate the
        highest score for each cell to create the matrix."""
        score_matrix = [[0 for j in range(cols)] for i in range(rows)]
        for i in range(rows):
            score_matrix[i][0] = -penalty * i
        for j in range(cols):
            score_matrix[0][j] = -penalty * j
        max_score = 0
        max_pos = None
        for i in range(1, rows):
            for j in range(1, cols):
                score = calc_score(score_matrix, i, j, matrix)
                if score > max_score:
                    max_score = score
                    max_pos = (i, j)
                score_matrix[i][j] = score
        assert max_pos is not None, 'the x, y position with the highest score was not found'
        return score_matrix

    def calc_score(score_matrix, i, j, matrix):
        """Calculates the score for each individual cell depending on the highest score from one of three directions:
        diagonal, up, and left. The scores are calculated by the following algorithms."""

        diag_score = score_matrix[i - 1][j - 1] + matrix[ord(seq1[i - 1]) - ord("A")][ord(seq2[j - 1]) - ord("A")]
        up_score = score_matrix[i - 1][j] - penalty
        left_score = score_matrix[i][j - 1] - penalty
        return max(diag_score, up_score, left_score)

    # Initialize the score matrix
    score_matrix = create_score_matrix(rows, cols, matrix, penalty)


    def traceback(score_matrix, matrix):
        """Traces back through the score_matrix comparing scores in each cell to possible scores from one of three
        directions: diagonal, up, and left. The highest score is assumed as correct, and if the current score matches
        the calculation for a direction that direction is followed to produce the alignment. Matches and gaps are
        assigned depending on direction, e.g. A gap in seq2 for a left movement and a gap in seq1 for an upwards movement.
        Two aligned sequences are constructed in this manner and returned."""
        aligned_seq1 = []
        aligned_seq2 = []
        n = len(seq1)
        m = len(seq2)
        i, j = n, m
        while i > 0 and j > 0:
            current = score_matrix[i][j]
            diagonal = score_matrix[i - 1][j - 1]
            up = score_matrix[i][j - 1]
            left = score_matrix[i - 1][j]
            if current == diagonal + matrix[ord(seq1[i - 1]) - ord("A")][ord(seq2[j - 1]) - ord("A")]:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif current == left - penalty:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            elif current == up - penalty:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j - 1])
                j -= 1

        while i > 0:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        while j > 0:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

        return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))

    # Initialize traceback function
    seq1_aligned, seq2_aligned = traceback(score_matrix, matrix)

    def final_align():
        def alignment_string(aligned_seq1, aligned_seq2):
            """This function adds in markers for matches, mismatches, and gaps and counts the numbers of each. A match is
            signified by a '|' symbol, a mismatch by ':', and a gap by '-', the numbers of each are counted and returned.
            The alignment score is also calculated by this function with exchange scores being added for matches and mismatches,
            and gap penalties added for the gap positions"""
            idents, gaps, mismatches = 0, 0, 0
            alignment_string = []
            score = 0
            for base1, base2 in zip(aligned_seq1, aligned_seq2):
                if base1 == base2:
                    alignment_string.append('|')
                    idents += 1
                    score += matrix[ord(base1) - ord("A")][ord(base2) - ord("A")]
                elif '-' in (base1, base2):
                    alignment_string.append(' ')
                    gaps += 1
                    score += -penalty
                else:
                    alignment_string.append(':')
                    mismatches += 1
                    score += matrix[ord(base1) - ord("A")][ord(base2) - ord("A")]
            return ''.join(alignment_string), idents, gaps, mismatches, score

        """This section of code prints the alignment in a similar way to a BLAST search output, it provides identities,
        gaps, shows the alignments with the identfiers seen in the previous function, and finally prints alignment score"""
        alignment_str, idents, gaps, mismatches, score = alignment_string(seq1_aligned, seq2_aligned)
        alength = len(seq1_aligned)
        print()
        print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents, alength, idents / alength, gaps,
                                                                             alength, gaps / alength))
        print()
        for i in range(0, alength, 60):
            seq1_slice = seq1_aligned[i:i + 60]
            print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
            print('             {0}'.format(alignment_str[i:i + 60]))
            seq2_slice = seq2_aligned[i:i + 60]
            print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
            print("Alignment score = %s" % score)
            print()



    alignment = final_align()
    return alignment, score_matrix


def do_local_alignment(sequences, matrix, penalty):
    """Do pairwise local alignment using DP."""
    """Assignment of variables for sequences, and the numbers of rows and columns for the score matrix"""
    seq1 = str(sequences[0].Sequence)
    seq2 = str(sequences[1].Sequence)
    rows = len(seq1) + 1
    cols = len(seq2) + 1


    def create_score_matrix(rows, cols, matrix, penalty):
        """This function creates the score matrix using the number of rows and columns assigned previously,
        as this is a local alignment algorithm the score matrix initializes with 0s in the first row and column.
        The calc_score function is included to fill in each cell independently"""
        score_matrix = [[0 for j in range(cols)] for i in range(rows)]
        max_score = 0
        max = None
        for i in range(1, rows):
            for j in range(1, cols):
                score = calc_score(score_matrix, i, j, matrix)
                if score > max_score:
                    max_score = score
                    max = (i, j)
                score_matrix[i][j] = score
        assert max is not None, 'The highest x, y position cannot be found'
        return score_matrix, max

    def calc_score(score_matrix, i, j, matrix):
        """This function calculates movement scores for each of the directions: diagonal, up, and left then returns
        the highest score. In this manner it is used in the create_score_matrix function to fill in each cell."""

        diag_score = score_matrix[i - 1][j - 1] + matrix[ord(seq1[i - 1]) - ord("A")][ord(seq2[j - 1]) - ord("A")]
        up_score = score_matrix[i - 1][j] - penalty
        left_score = score_matrix[i][j - 1] - penalty
        return max(0, diag_score, up_score, left_score)



    # Initialize the score matrix
    score_matrix, start_pos = create_score_matrix(rows, cols, matrix, penalty)

    def traceback(score_matrix, matrix):
        """The traceback function follows the alignment path, as this is a local alignment the traceback starts at the
        highest scoring cell in the highest scoring cell. Two lists representing the aligned sequences are appended
        with matches or gaps as needed."""
        aligned_seq1 = []
        aligned_seq2 = []
        n = len(seq1)
        m = len(seq2)
        i, j = start_pos
        score = 0
        while i > 0 and j > 0:
            current = score_matrix[i][j]
            diagonal = score_matrix[i - 1][j - 1]
            up = score_matrix[i][j - 1]
            left = score_matrix[i - 1][j]
            if current == diagonal + matrix[ord(seq1[i - 1]) - ord("A")][ord(seq2[j - 1]) - ord("A")]:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif current == left - penalty:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            elif current == up - penalty:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j - 1])
                j -= 1

        while i > 0:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        while j > 0:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

        return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))




        # Initialize traceback function

    seq1_aligned, seq2_aligned = traceback(score_matrix, matrix)
    def final_align():
        def alignment_string(aligned_seq1, aligned_seq2):
            """This function creates a string between the two aligned sequences that provides a visual identifier
            for matches, mismatches, and gaps which are '|', ':', and ' ' respectively. It also calculates the numbers
            of each of these for later use and calculates the alignment score."""
            idents, gaps, mismatches = 0, 0, 0
            alignment_string = []
            score = 0
            for base1, base2 in zip(aligned_seq1, aligned_seq2):
                if base1 == base2:
                    alignment_string.append('|')
                    idents += 1
                    score += matrix[ord(base1) - ord("A")][ord(base2) - ord("A")]
                elif '-' in (base1, base2):
                    alignment_string.append(' ')
                    gaps += 1
                    score += -penalty
                else:
                    alignment_string.append(':')
                    mismatches += 1
                    score += matrix[ord(base1) - ord("A")][ord(base2) - ord("A")]
            return ''.join(alignment_string), idents, gaps, mismatches, score


        """Prints finished alignment in a similar manner to a BLAST output"""
        alignment_str, idents, gaps, mismatches, score = alignment_string(seq1_aligned, seq2_aligned)
        alength = len(seq1_aligned)
        print()
        print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents, alength, idents / alength, gaps,
                                                                             alength, gaps / alength))
        print()
        for i in range(0, alength, 60):
            seq1_slice = seq1_aligned[i:i + 60]
            print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
            print('             {0}'.format(alignment_str[i:i + 60]))
            seq2_slice = seq2_aligned[i:i + 60]
            print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
            print("Alignment score = %s" % score)
            print()





    alignment = final_align()
    return alignment, score_matrix



def do_semiglobal_alignment(sequences, matrix, penalty):
    """Do pairwise semi-global alignment using DP."""
    """Set inital variables for the two sequences and for the number of rows and columns needed for the matrix"""
    seq1 = str(sequences[0].Sequence)
    seq2 = str(sequences[1].Sequence)
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    def create_score_matrix(rows, cols, matrix, penalty):
        """Creates the score matrix using the number of rows and columns assigned previously,
        as this is a semi-global alignment initial gaps are free and as such the matrix initializes
        with zeroes in the first column and row. Again the calc_score function is used to calculate the
        highest score for each cell to create the matrix."""
        score_matrix = [[0 for j in range(cols)] for i in range(rows)]
        max_score = 0
        max = None
        for i in range(1, rows):
            for j in range(1, cols):
                score = calc_score(score_matrix, i, j, matrix)
                if score > max_score:
                    max_score = score
                    max = (i, j)
                score_matrix[i][j] = score
        assert max is not None, 'The highest x, y position cannot be found'
        return score_matrix, max

    def calc_score(score_matrix, i, j, matrix):
        """Calculates the score for each individual cell depending on the highest score from one of three directions:
        diagonal, up, and left. The scores are calculated by the following algorithms."""

        diag_score = score_matrix[i - 1][j - 1] + matrix[ord(seq1[i - 1]) - ord("A")][ord(seq2[j - 1]) - ord("A")]
        up_score = score_matrix[i - 1][j] - penalty
        left_score = score_matrix[i][j - 1] - penalty
        return max(diag_score, up_score, left_score)

    # Initialize the score matrix


    score_matrix, start_pos = create_score_matrix(rows, cols, matrix, penalty)


    def traceback(score_matrix, matrix):
        """Traces back through the score_matrix comparing scores in each cell to possible scores from one of three
        directions: diagonal, up, and left. The highest score is assumed as correct, and if the current score matches
        the calculation for a direction that direction is followed to produce the alignment. Matches and gaps are
        assigned depending on direction, e.g. A gap in seq2 for a left movement and a gap in seq1 for an upwards movement.
        Two aligned sequences are constructed in this manner and returned."""
        aligned_seq1 = []
        aligned_seq2 = []
        n = len(seq1)
        m = len(seq2)
        i, j = n, m
        while i > 0 and j > 0:
            current = score_matrix[i][j]
            diagonal = score_matrix[i - 1][j - 1]
            up = score_matrix[i][j - 1]
            left = score_matrix[i - 1][j]
            if current == diagonal + matrix[ord(seq1[i - 1]) - ord("A")][ord(seq2[j - 1]) - ord("A")]:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif current == left - penalty:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            elif current == up - penalty:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j - 1])
                j -= 1

        while i > 0:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        while j > 0:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

        return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))




        # Initialize traceback function

    seq1_aligned, seq2_aligned = traceback(score_matrix, matrix)


    def final_align():
        def alignment_string(aligned_seq1, aligned_seq2):
            """This function assigns visual identifiers to match, mismatch, and gap placements which are '|', ':', and '-'
            respectively """
            idents, gaps, mismatches = 0, 0, 0
            alignment_string = []
            score = 0
            for base1, base2 in zip(aligned_seq1, aligned_seq2):
                if base1 == base2:
                    alignment_string.append('|')
                    idents += 1
                    score += matrix[ord(base1) - ord("A")][ord(base2) - ord("A")]
                elif '-' in (base1, base2):
                    alignment_string.append(' ')
                    gaps += 1
                    score += -penalty
                else:
                    alignment_string.append(':')
                    mismatches += 1
                    score += matrix[ord(base1) - ord("A")][ord(base2) - ord("A")]
            if '-' in (aligned_seq1[len(aligned_seq1) - 1], aligned_seq2[len(aligned_seq2) - 1]):
                    score += penalty
            return ''.join(alignment_string), idents, gaps, mismatches, score

        alignment_str, idents, gaps, mismatches, score = alignment_string(seq1_aligned, seq2_aligned)
        alength = len(seq1_aligned)
        print()
        print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents, alength, idents / alength, gaps,
                                                                             alength, gaps / alength))
        print()
        for i in range(0, alength, 60):
            seq1_slice = seq1_aligned[i:i + 60]
            print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
            print('             {0}'.format(alignment_str[i:i + 60]))
            seq2_slice = seq2_aligned[i:i + 60]
            print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
            print("Alignment score = %s" % score)
            print()



    alignment = final_align()
    return alignment, score_matrix


def print_matrix_to_file(score_matrix, fileName):
    """Write a matrix into file.

    matrix: a list of list in Python, storing an alignment or a score
    matrix.
    fileName: str, a file name (with a path) to store the matrix.
    It is not recommended to tinker with this function.
    """
    with open(fileName, 'w') as f:
        for row in score_matrix:
            print(''.join(map(str, row)), file=f)


def print_alignment_to_file(alignment, fileName):
    """Writes alignment to specified output file using the created alignment and specified output file"""

    with open(fileName, 'w') as f:
        print(alignment)


def print_matrix_on_screen(score_matrix, width=5):
    """Print a matrix on the screen.

    matrix: a list of list in Python, storing an alignment or a score
    matrix.
    width: that of the space one cell occupies.
    This will facilitate your testing.
    """
    for row in score_matrix:
        print(''.join(['{0:>{w}}'.format(item, w=width) for item in row]))


def main():
    """Main function.

    Please change it accordingly to make the program work.
    """
    # get command line options
    global sequences
    args = get_args()

    # set substitution matrix:
    if args.exchange_matrix == "pam250":
        exchangeMatrix = pam250
    elif args.exchange_matrix == "blosum62":
        exchangeMatrix = blosum62
    else:
        exchangeMatrix = identity

    print("Exchange matrix is %s" % args.exchange_matrix.upper())
    print("Gap penalty = %s" % args.gap_penalty)

    # read sequences from fasta file, and catch error reading file
    try:
        sequences = readSequences(open(args.fasta))
    except OSError as e:
        print("ERROR: cannot open or read fasta input file:", e.filename)

    for seq in sequences:
        print(seq)
        print("Sequence length = %s" % len(seq.Sequence))

    # call alignment routine(s):
    if args.align_global:
        alignment, score_matrix = do_global_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    elif args.align_local:
        alignment, score_matrix = do_local_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    elif args.align_semiglobal:
        alignment, score_matrix = do_semiglobal_alignment(
                sequences, exchangeMatrix, args.gap_penalty)
    else:
        sys.exit("BUG! this should not happen.")

    if args.alignment:
        print_alignment_to_file(alignment, args.alignment)
    if args.score_matrix:
        print_matrix_to_file(score_matrix, args.score_matrix)
    if args.print_on_screen:
        print_matrix_on_screen(score_matrix, width=5)


if __name__ == "__main__":
    main()

# last line