import numpy as np
from xlwings import xrange

global match
global mismatch
global gap_open
global gap_extension
# global M
# global Ix
# global Iy
Infinity = float('inf')

match = 1
mismatch = -1
gap_open = -10
gap_extension = -1


def split(sequence):
    return [char for char in sequence]


def make_matrix(row, cals):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0] * cals for i in range(row)]


class Affine_gap:
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2

    def matchchar(self, seq1, seq2, i, j):
        if seq2[i - 1] == seq1[j - 1]:
            return match
        else:
            return mismatch

    def alingment(self):
        dim_seq1 = len(self.seq1) + 1
        split_seq1 = split(self.seq1)
        dime_seq2 = len(self.seq2) + 1
        split_seq2 = split(self.seq2)
        M = make_matrix(dime_seq2, dim_seq1)
        Ix = make_matrix(dime_seq2, dim_seq1)
        Iy = make_matrix(dime_seq2, dim_seq1)
        M[0][0] = 0
        Ix[0][0] = gap_open - gap_extension
        Iy[0][0] = gap_open - gap_extension

        for i in range(1, dime_seq2):
            M[i][0] = -Infinity
            Ix[i][0] = gap_open + ((i - 1) * gap_extension)
            Iy[i][0] = -Infinity

        for j in range(1, dim_seq1):
            M[0][j] = -Infinity
            Ix[0][j] = -Infinity
            Iy[0][j] = gap_open + ((j - 1) * gap_extension)

        for i in range(1, dime_seq2):
            for j in range(1, dim_seq1):
                M[i][j] = self.matchchar(self.seq1, self.seq2, i, j) + max(
                    M[i - 1][j - 1],
                    Ix[i - 1][j - 1],
                    Iy[i - 1][j - 1]
                )

                Ix[i][j] = max(
                    gap_open + M[i - 1][j],
                    gap_extension + Ix[i - 1][j],
                )

                Iy[i][j] = max(
                    gap_open + M[i][j - 1],
                    gap_extension + Iy[i][j - 1]
                )

        Ix_score = (Ix[i][j])
        Iy_score = (Iy[i][j])
        M_score = (M[i][j])
        # print(M)
        # print(Ix)
        # print(Iy)
        optimal_score = max(Ix_score, Iy_score, M_score)
        seq1 = split_seq1
        seq2 = split_seq2
        GAFirst = ""
        GASecond = ""
        if M[i][j] > Ix[i][j] and M[i][j] > Ix[i][j]:
            current_matrix = 'm'
            optimalScore = M[i][j]

        elif Iy[i][j] > Ix[i][j] and Iy[i][j] > M[i][j]:
            current_matrix = 'y'
            optimalScore = Iy[i][j]

        else:
            current_matrix = 'x'
            optimalScore = Ix[i][j]

        while i > 0 or j > 0:
            if current_matrix == 'm':
                GAFirst += seq1[j - 1]
                GASecond += seq2[i - 1]
                if seq2[i - 1] == seq1[j - 1]:
                    penalty = match
                else:
                    penalty = mismatch

                if (M[i - 1][j - 1] + penalty) == M[i][j]:
                    i -= 1
                    j -= 1
                    current_matrix = 'm'
                elif (Ix[i - 1][j - 1] + penalty) == M[i][j]:
                    i -= 1
                    j -= 1
                    current_matrix = 'x'
                elif (Iy[i - 1][j - 1] + penalty) == M[i][j]:
                    i -= 1
                    j -= 1
                    current_matrix = 'y'

            elif current_matrix == 'x':
                GAFirst += "-"
                GASecond += seq2[i - 1]

                if (Ix[i - 1][j] + gap_extension) == Ix[i][j]:
                    i -= 1
                    current_matrix = 'x'
                elif (M[i - 1][j] + gap_open) == Ix[i][j]:
                    i -= 1
                    current_matrix = 'm'

            elif current_matrix == 'y':
                GAFirst += seq1[j - 1]
                GASecond += "-"
                if (Iy[i][j - 1] + gap_extension) == Iy[i][j]:
                    j -= 1
                    current_matrix = 'y'
                elif (M[i][j - 1] + gap_open) == Iy[i][j]:
                    j -= 1
                    current_matrix = 'm'
        GAFirst = GAFirst[::-1]
        GASecond = GASecond[::-1]
        return GAFirst, GASecond, optimal_score


input_seq = []
for i in xrange(4):
    try:
        line = input()
    except EOFError:
        break
    input_seq.append(line)

seq_0_1 = Affine_gap(input_seq[0], input_seq[1])
GA_first_seq_0_1, GA_second_seq_0_1, num_seq_0_1 = seq_0_1.alingment()
# print(num_seq_0_1)

seq_0_2 = Affine_gap(input_seq[0], input_seq[2])
GA_first_seq_0_2, GA_second_seq_0_2, num_seq_0_2 = seq_0_2.alingment()
# print(num_seq_0_2)

seq_0_3 = Affine_gap(input_seq[0], input_seq[3])
GA_first_seq_0_3, GA_second_seq_0_3, num_seq_0_3 = seq_0_3.alingment()
# print(num_seq_0_3)

seq_1_0 = Affine_gap(input_seq[1], input_seq[0])
GA_first_seq_1_0, GA_second_seq_1_0, num_seq_1_0 = seq_1_0.alingment()
# print(num_seq_1_0)

seq_1_2 = Affine_gap(input_seq[1], input_seq[2])
GA_first_seq_1_2, GA_second_seq_1_2, num_seq_1_2 = seq_1_2.alingment()
# print(num_seq_1_2)

seq_1_3 = Affine_gap(input_seq[1], input_seq[3])
GA_first_seq_1_3, GA_second_seq_1_3, num_seq_1_3 = seq_1_3.alingment()
# print(num_seq_1_3)

seq_2_0 = Affine_gap(input_seq[2], input_seq[0])
GA_first_seq_2_0, GA_second_seq_2_0, num_seq_2_0 = seq_2_0.alingment()
# print(num_seq_2_0)

seq_2_1 = Affine_gap(input_seq[2], input_seq[1])
GA_first_seq_2_1, GA_second_seq_2_1, num_seq_2_1 = seq_2_1.alingment()
# print(num_seq_2_1)

seq_2_3 = Affine_gap(input_seq[2], input_seq[3])
GA_first_seq_2_3, GA_second_seq_2_3, num_seq_2_3 = seq_2_3.alingment()
# print(num_seq_2_3)

seq_3_0 = Affine_gap(input_seq[3], input_seq[0])
GA_first_seq_3_0, GA_second_seq_3_0, num_seq_3_0 = seq_3_0.alingment()
# print(num_seq_3_0)

seq_3_1 = Affine_gap(input_seq[3], input_seq[1])
GA_first_seq_3_1, GA_second_seq_3_1, num_seq_3_1 = seq_3_1.alingment()
# print(num_seq_3_1)

seq_3_2 = Affine_gap(input_seq[3], input_seq[2])
GA_first_seq_3_2, GA_second_seq_3_2, num_seq_3_2 = seq_3_2.alingment()
# print(num_seq_3_2)

sum_s0 = num_seq_0_1 + num_seq_0_2 + num_seq_0_3
sum_s1 = num_seq_1_0 + num_seq_1_2 + num_seq_1_3
sum_s2 = num_seq_2_0 + num_seq_2_1 + num_seq_2_3
sum_s3 = num_seq_3_0 + num_seq_3_1 + num_seq_3_2

# print(sum_s0)
# print(sum_s1)
# print(sum_s2)
# print(sum_s3)

sum_list = [sum_s0, sum_s1, sum_s2, sum_s3]
max_val = max(sum_s0, sum_s1, sum_s2, sum_s3)
# print(max_val)
if sum_list[0] == max_val:
    # print("The Sc is first sequence: ")
    res_seq_0_2 = Affine_gap(GA_second_seq_0_1, GA_second_seq_0_2)
    GA_second_seq_0_1, GA_second_seq_0_2, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_second_seq_0_1, GA_second_seq_0_3)
    GA_second_seq_0_1, GA_second_seq_0_3, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_second_seq_0_2, GA_second_seq_0_3)
    GA_second_seq_0_2, GA_second_seq_0_3, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_first_seq_0_3, GA_second_seq_0_3)
    GA_first_seq_0_3, GA_second_seq_0_3, num_seq_1_0 = res_seq_0_2.alingment()

    print(GA_second_seq_0_1)
    print(GA_first_seq_0_3)
    print(GA_second_seq_0_2)
    print(GA_second_seq_0_3)

if sum_list[1] == max_val:
    # print("The Sc is second sequence: ")

    res_seq_0_2 = Affine_gap(GA_second_seq_1_0, GA_second_seq_1_2)
    GA_second_seq_1_0, GA_second_seq_1_2, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_second_seq_1_0, GA_second_seq_1_3)
    GA_second_seq_1_0, GA_second_seq_1_3, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_second_seq_1_2, GA_second_seq_1_3)
    GA_second_seq_1_2, GA_second_seq_1_3, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_first_seq_1_0, GA_second_seq_1_3)
    GA_first_seq_1_0, GA_second_seq_0_3, num_seq_1_0 = res_seq_0_2.alingment()

    print(GA_second_seq_1_0)
    print(GA_first_seq_1_0)
    print(GA_second_seq_1_2)
    print(GA_second_seq_1_3)

if sum_list[2] == max_val:
    # print("The Sc is third sequence: ")

    res_seq_0_2 = Affine_gap(GA_second_seq_2_0, GA_second_seq_2_1)
    GA_second_seq_2_0, GA_second_seq_2_1, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_second_seq_2_0, GA_second_seq_2_3)
    GA_second_seq_2_0, GA_second_seq_2_3, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_second_seq_2_1, GA_second_seq_2_3)
    GA_second_seq_2_1, GA_second_seq_2_3, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_first_seq_2_0, GA_second_seq_2_3)
    GA_first_seq_2_0, GA_second_seq_2_3, num_seq_1_0 = res_seq_0_2.alingment()

    print(GA_second_seq_2_0)
    print(GA_second_seq_2_1)
    print(GA_first_seq_2_0)
    print(GA_second_seq_2_3)

if sum_list[3] == max_val:
    # print("The Sc is fourth sequence: ")

    res_seq_0_2 = Affine_gap(GA_second_seq_3_0, GA_second_seq_3_1)
    GA_second_seq_3_0, GA_second_seq_3_1, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_second_seq_3_0, GA_second_seq_3_2)
    GA_second_seq_3_0, GA_second_seq_3_2, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_second_seq_3_1, GA_second_seq_3_2)
    GA_second_seq_3_1, GA_second_seq_3_2, num_seq_1_0 = res_seq_0_2.alingment()
    res_seq_0_2 = Affine_gap(GA_first_seq_3_0, GA_second_seq_3_2)
    GA_first_seq_3_0, GA_second_seq_3_2, num_seq_1_0 = res_seq_0_2.alingment()

    print(GA_second_seq_3_0)
    print(GA_second_seq_3_1)
    print(GA_second_seq_3_2)
    print(GA_first_seq_3_0)
