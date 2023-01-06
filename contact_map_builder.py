import sys
import numpy as np

from collections import defaultdict


def read_msa(file_name):
    msa = []
    start_idx = None
    with open(file_name) as f:
        for l in f.readlines():
            # once need to find the place in the line where the sequence starts
            # if there is no spaces in the line then format is wrong
            if start_idx is None:
                space_idx = l.rfind(' ')
                if space_idx == -1:
                    raise TypeError('File format is wrong. Should be "Selex".')
                start_idx = space_idx + 1

            msa.append(list(l[start_idx:].strip('\n')))
    return np.array(msa)


def save_cm_ft_comar_format(cm, output_file_path):
    with open(output_file_path, 'w') as f:
        f.write(f'{cm.shape[0]}\n')
        f.write('\n'.join([str(x) for x in cm.reshape(-1)]))


def calculate_residue_probabilities(msa):
    M, N = msa.shape
    p = []
    for i in range(N):
        p_i = defaultdict(float)
        for l in msa[:, i]:
            p_i[l] += 1 / M
        p.append(p_i)
    return p


def calculate_pair_probabilities(msa):
    M, N = msa.shape
    p = []
    for i in range(N):
        p_i = []
        for j in range(N):
            p_ij = defaultdict(float)
            for l_i, l_j in zip(msa[:, i], msa[:, j]):
                p_ij[l_i + l_j] += 1 / M
            p_i.append(p_ij)
        p.append(p_i)
    return p


def calculate_MI_matrix(msa, p, p_pair):
    M, N = msa.shape
    mi_matrix = np.empty((N, N))
    msa_abc = set(msa.ravel()) - {'.'}

    for i in range(N):
        for j in range(N):
            sum_ab = 0
            for a in msa_abc:
                for b in msa_abc:
                    p_a = p[i][a]
                    p_b = p[j][b]
                    p_ab = p_pair[i][j][a + b]

                    if p_a != 0 and p_b != 0 and p_ab != 0:
                        sum_ab += p_ab * np.log(p_ab / (p_a * p_b))
            mi_matrix[i, j] = sum_ab
    return mi_matrix


def do_MI_correction(mi_matrix):
    N = mi_matrix.shape[0]
    mi_matrix_correction = np.empty((N, N))

    for i in range(N):
        for j in range(N):
            correction = np.sum([(mi_matrix[k, j] + mi_matrix[i, k]) for k in range(N)]) / N
            mi_matrix_correction[i, j] = mi_matrix[i, j] - correction
    
    return mi_matrix_correction


def main(msa_file_path, output_cm_file_path):
    # assume that this threshold is good for all msa
    threshold = 0.01

    # read msa from Selex format file
    msa = read_msa(msa_file_path)

    # calculate probabilities
    p = calculate_residue_probabilities(msa)
    p_pair = calculate_pair_probabilities(msa)

    # calculate mutual information matrix and correct it
    mi_matrix = calculate_MI_matrix(
        msa=msa, p=p, p_pair=p_pair
    )
    mi_matrix_correction = do_MI_correction(mi_matrix)

    # built contact map by using hardcoded threshold
    cm = np.int0(mi_matrix_correction > threshold)

    # save contact map with readable for FT-comar format
    save_cm_ft_comar_format(cm, output_cm_file_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])