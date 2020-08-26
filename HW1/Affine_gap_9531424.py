import numpy as np

Infinity = float('inf')


def split(sequence):
    return [char for char in sequence]


def make_matrix(row, cals):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0] * cals for i in range(row)]


class Affine_gap:
    def __init__(self, seq1, seq2, match, mismatch, gap_open, gap_extension):
        self.seq1 = seq1
        self.seq2 = seq2
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extension = gap_extension
        self.alingment()

    def matchchar(self, seq1, seq2, i, j):
        if seq2[i - 1] == seq1[j - 1]:
            return self.match
        else:
            return self.mismatch

    def alingment(self):
        dim_seq1 = len(self.seq1) + 1
        split_seq1 = split(self.seq1)
        dime_seq2 = len(self.seq2) + 1
        split_seq2 = split(self.seq2)
        M = make_matrix(dime_seq2, dim_seq1)
        Ix = make_matrix(dime_seq2, dim_seq1)
        Iy = make_matrix(dime_seq2, dim_seq1)
        M[0][0] = 0
        Ix[0][0] = self.gap_open - self.gap_extension
        Iy[0][0] = self.gap_open - self.gap_extension

        for i in range(1, dime_seq2):
            M[i][0] = -Infinity
            Ix[i][0] = self.gap_open + ((i - 1) * self.gap_extension)
            Iy[i][0] = -Infinity

        for j in range(1, dim_seq1):
            M[0][j] = -Infinity
            Ix[0][j] = -Infinity
            Iy[0][j] = self.gap_open + ((j - 1) * self.gap_extension)

        for i in range(1, dime_seq2):
            for j in range(1, dim_seq1):
                M[i][j] = self.matchchar(self.seq1, self.seq2, i, j) + max(
                    M[i - 1][j - 1],
                    Ix[i - 1][j - 1],
                    Iy[i - 1][j - 1]
                )

                Ix[i][j] = max(
                    self.gap_open + M[i - 1][j],
                    self.gap_extension + Ix[i - 1][j],
                )

                Iy[i][j] = max(
                    self.gap_open + M[i][j - 1],
                    self.gap_extension + Iy[i][j - 1]
                )
        Ix_score = (Ix[i][j])
        Iy_score = (Iy[i][j])
        M_score = (M[i][j])
        optimal_score = max(Ix_score, Iy_score, M_score)
        print(optimal_score)
        self.traceback(M, Ix, Iy, self.gap_open, self.gap_extension, split_seq1, split_seq2, i, j)

    def traceback(self, M, Ix, Iy, gap_open, gap_extension, seq1, seq2, row, col):
        GAFirst = ""
        GASecond = ""
        i = row
        j = col
        if M[i][j] > Ix[i][j] and M[i][j] > Ix[i][j]:
            current_matrix = 'm'
            optimalScore = M[i][j]
            print("Optimal Score for M is  = " + str(optimalScore))

        elif Iy[i][j] > Ix[i][j] and Iy[i][j] > M[i][j]:
            current_matrix = 'y'
            optimalScore = Iy[i][j]
            print("Optimal Score for Iy is = " + str(optimalScore))

        else:
            current_matrix = 'x'
            optimalScore = Ix[i][j]
            print("Optimal Score for Ix is = " + str(optimalScore))

        while i > 0 or j > 0:
            if current_matrix == 'm':
                GAFirst += seq1[j - 1]
                GASecond += seq2[i - 1]
                print(GASecond)
                if seq2[i - 1] == seq1[j - 1]:
                    penalty = self.match
                else:
                    penalty = self.mismatch

                if (Iy[i - 1][j - 1] + penalty) == M[i][j]:
                    print("section 1")
                    print(M[i - 1][j - 1])
                    print(penalty)
                    print(M[i][j])
                    print("End section 1")

                    i -= 1
                    j -= 1
                    current_matrix = 'y'
                elif (Ix[i - 1][j - 1] + penalty) == M[i][j]:
                    print("section 2")
                    print(Ix[i - 1][j - 1])
                    print(penalty)
                    print(M[i][j])
                    print("End section 2")
                    i -= 1
                    j -= 1
                    current_matrix = 'x'
                elif (M[i - 1][j - 1] + penalty) == M[i][j]:
                    print("section 3")
                    print(Iy[i - 1][j - 1])
                    print(penalty)
                    print(M[i][j])
                    print("End section 3")
                    i -= 1
                    j -= 1
                    current_matrix = 'm'

            elif current_matrix == 'x':
                GAFirst += "-"
                GASecond += seq2[i - 1]

                if (Ix[i - 1][j] + gap_extension) == Ix[i][j]:
                    print("Sub section 1 for X: ")
                    print(Ix[i - 1][j])
                    print(gap_extension)
                    print(Ix[i][j])
                    i -= 1
                    current_matrix = 'x'
                elif (M[i - 1][j] + gap_open) == Ix[i][j]:
                    print("Sub section 2 for X: ")
                    print(M[i - 1][j])
                    print(gap_open)
                    print(Ix[i][j])
                    i -= 1
                    current_matrix = 'm'

            elif current_matrix == 'y':
                GAFirst += seq1[j - 1]
                GASecond += "-"

                if (Iy[i][j - 1] + gap_extension) == Iy[i][j]:
                    print("Sub section 1 for Y: ")
                    print(Iy[i][j - 1])
                    print(gap_extension)
                    print(Iy[i][j])
                    j -= 1
                    current_matrix = 'y'
                elif (M[i][j - 1] + gap_open) == Iy[i][j]:
                    print("Sub section 2 for Y: ")
                    print(M[i][j - 1])
                    print(gap_open)
                    print(Iy[i][j])
                    j -= 1
                    current_matrix = 'm'
        GAFirst = GAFirst[::-1]
        GASecond = GASecond[::-1]

        print("Optimal Alignment:")
        print(GAFirst)
        print(GASecond)


input_seq = ['ACACT', 'AAT']

model = Affine_gap(input_seq[0], input_seq[1], 1, -1, -4, -1)
