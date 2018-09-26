"""Script for computing sequence alignments using Needleman-Wunsch with
   affine gap penalties.
Arguments:
    f - FASTA file with sequences in FASTA format.
    s - JSON with the score matrix for alignment.
    d - The gap opening penalty for the alignment.
    e - The gap extension penalty for the alignment.

Outputs:
    Prints alignment to console.

Example Usage:
    python 2c.py -f sequences.fasta -s score_matrix.json -d 430 -e 30
"""

import argparse
import json
import numpy as np

'''Computes the actual string alignments given the traceback matrix.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    t: the traceback matrix, which stores values that point to which
       prior matrix was used to reach a given location in each of the
       3 matrices.
    start: value indicating the starting matrix (that had the optimal value)
Returns:
    a_x: the string for the alignment of x's sequence
    a_y: the string for the alignment of y's sequence
'''


def traceback(x, y, t, start):
    # Initialize a_x and a_y
    a_x = ''
    a_y = ''

    # Starting from the bottom right corner of the matrix to work backwards
    r = len(x)
    c = len(y)

    # Set current value to predetermined starting point
    curr = start

    # Until reaching top left corner, add next residue or gap to a_x or a_y
    while r != 0 or c != 0:
        if curr == 0:
            a_x += x[r-1]
            a_y += y[c-1]
            r -= 1
            c -= 1
        elif curr == 1:
            a_x += x[r-1]
            a_y += '-'
            r -= 1
        else:
            a_x += '-'
            a_y += y[c-1]
            c -= 1

        # Switch between m, i_x, and i_y if necessary before next iteration
        prev = t[r][c][curr]
        curr = prev

    # Reverse strings (we were appending to the end not beginning) and return
    return a_x[::-1], a_y[::-1]


'''Computes the score and alignment of two strings using an affine gap penalty.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    s: the score matrix
    d: the gap opening penalty
    e: the gap extension penalty
Returns:
    score: the score of the optimal sequence alignment
    a_x: the aligned first string
    a_y: the aligned second string
The latter two are computed using the above traceback method.
'''


def affine_sequence_alignment(x, y, s, d, e):

    # Initialize necessary matrices
    m = np.zeros([len(x)+1, len(y)+1])
    ''' Recurrence matrix'''

    i_x = np.zeros([len(x)+1, len(y)+1])
    ''' Recurrence matrix'''

    i_y = np.zeros([len(x)+1, len(y)+1])
    ''' Recurrence matrix'''

    # Extra dimension
    t = np.zeros([len(x)+1, len(y)+1, 3])
    ''' Traceback matrix, redefine/use as necessary. '''


    for i in range(1, len(x)+1):
        m[i][0] = float("-inf")
        i_x[i][0] = i_x[i-1][0] - e
        i_y[i][0] = float("-inf")

    for j in range(1, len(y)+1):
        m[0][j] = float("-inf")
        i_x[0][j] = float("-inf")
        i_y[0][j] = i_y[0][j-1]-e

    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):

            m_val = max(m[i-1][j-1]+s[x[i-1]][y[j-1]],
                        i_x[i-1][j-1]+s[x[i-1]][y[j-1]],
                        i_y[i-1][j-1]+s[x[i-1]][y[j-1]])
            m[i][j] = m_val
            if m_val == m[i-1][j-1]+s[x[i-1]][y[j-1]]:
                t[i][j][0] = 0
            elif m_val == i_x[i-1][j-1]+s[x[i-1]][y[j-1]]:
                t[i][j][0] = 1
            else:
                t[i][j][0] = 2

            ix_val = max(m[i-1, j]-d,
                         i_x[i-1][j]-e)
            i_x[i][j] = ix_val
            if ix_val == m[i-1, j]-d:
                t[i][j][1] = 0
            else:
                t[i][j][1] = 1

            iy_val = max(m[i][j-1]-d,
                         i_y[i][j-1]-e)
            i_y[i][j] = iy_val
            if iy_val == m[i][j-1]-d:
                t[i][j][2] = 0
            else:
                t[i][j][2] = 2

    score = max(m[len(x)][len(y)],
                i_x[len(x)][len(y)],
                i_y[len(x)][len(y)])

    if score == m[len(x)][len(y)]:
        start = 0
    elif score == i_x[len(x)][len(y)]:
        start = 1
    else:
        start = 2
    ''' Indicator of the starting point (which matrix has optimal value).'''
    a_x, a_y = traceback(x, y, t, start)
    return score, (a_x, a_y)


'''Prints two aligned sequences formatted for convenient inspection.
Arguments:
    a_x: the first sequence aligned
    a_y: the second sequence aligned
Outputs:
    Prints aligned sequences (80 characters per line) to console
'''


def print_alignment(a_x, a_y):
    assert len(a_x) == len(a_y), "Sequence alignment lengths must be the same."
    for i in range(1 + (len(a_x) / 80)):
        start = i * 80
        end = (i + 1) * 80
        print a_x[start:end]
        print a_y[start:end]
        print


def main():
    parser = argparse.ArgumentParser(
        description='Calculate sequence alignments for two sequences with an affine gap penalty.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-s', action="store", dest="s", type=str, required=True)
    parser.add_argument('-d', action="store", dest="d", type=float, required=True)
    parser.add_argument('-e', action="store", dest="e", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    score_matrix_file = args.s
    d = args.d
    e = args.e

    with open(fasta_file) as f:
        _, x, _, y = [line.strip() for line in f.readlines()]
    with open(score_matrix_file) as f:
        s = json.loads(f.readlines()[0])

    score, (a_x, a_y) = affine_sequence_alignment(x, y, s, d, e)
    print "Alignment:"
    print_alignment(a_x, a_y)
    print "Score: " + str(score)


if __name__ == "__main__":
    main()
