from collections import Counter
from itertools import combinations_with_replacement as cwr
import numpy as np


# initial error rates
error_rates = {
    "A": {"A": 0.99,   "C": 0.0025, "G": 0.0025, "T": 0.0025, "-": 0.0025},
    "C": {"A": 0.0025, "C": 0.99,   "G": 0.0025, "T": 0.0025, "-": 0.0025},
    "G": {"A": 0.0025, "C": 0.0025, "G": 0.99,   "T": 0.0025, "-": 0.0025},
    "T": {"A": 0.0025, "C": 0.0025, "G": 0.0025, "T": 0.99,   "-": 0.0025},
    "-": {"A": 0,      "C": 0,      "G": 0,      "T": 0,      "-": 1},
}


# bases order in matrix ACGT-
bases_index = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "-": 4
}


alignment = np.array([
    ["A", "T", "A", "A", "C", "G", "T"],
    ["A", "T", "A", "A", "C", "G", "T"],
    ["A", "T", "A", "A", "C", "G", "T"],
    ["A", "T", "A", "A", "C", "G", "T"],
    ["A", "T", "A", "A", "C", "G", "T"],
    ["A", "T", "A", "A", "A", "G", "T"],
    ["A", "C", "A", "A", "A", "G", "T"],
    ["A", "G", "A", "A", "A", "G", "T"],
    ["A", "G", "A", "A", "A", "G", "T"],
    ["A", "G", "A", "A", "A", "G", "T"],
    ["A", "G", "A", "A", "C", "G", "T"],
    ["A", "G", "A", "A", "C", "G", "T"],
])


def all_genotype(ploidy):
    """
    int ploidy
    return all possible genotypes, completely determined by ploidy
    """
    return ["".join(comb) for comb in cwr("ACGT-", ploidy)]


def prob_t_N(genotype, base):
    """
    str genotype
    str base
    return P(base in genotype)

    diploid:    P_AA(A) = 1, P_AA(C) = 0, P_AC(A) = 0.5, P_AC(C) = 0.5
    triploid:   P_AAA(A) = 1, P_AAC(A) = 2/3
    tetraploid: P_AAAA(A) = 1, P_AACC(A) = 0.5, P_ACGT(A) = 0.25
    """
    cnter = Counter(genotype)
    return cnter.get(base, 0) * 1/len(genotype)


def likelihood_genotype(genotype, bases_all_reads):
    """
    str genotype
    iterable-obj bases_all_reads, list or np.array
    return P(data|genotype) == likelihood

    bases_all_reads: all bases in read at a given site
    """
    likelihood = 1
    for observed_base in bases_all_reads:
        p = 0
        for base in "ACGT-":
            l = prob_t_N(genotype, base) * error_rates[base][observed_base]
            p += l
        likelihood *= p

    return likelihood


def prob_N_data(bases_all_reads):
    """
    prod: prod of P(rih=N, nih|s) == P(rih=N,nih|s) * prod of P(njh|s)
    m_small: P(N|geno) of bases_all_reads at a given site
    m_big: P()
    """
    m_small = np.zeros((genotype_num, 5))
    m_big = np.zeros((K, 5))

    for i, observed_base in enumerate(bases_all_reads):
        rest_bases = bases_all_reads[:i] + bases_all_reads[i+1:]
        for b, possible_base in enumerate("ACGT-"):
            for g, geno in enumerate(all_geno):
                prod_i = (prob_t_N(geno, possible_base) * 
                          error_rates[possible_base][observed_base])
                prod_j = 1
                for j in rest_bases:
                    p_j = 0
                    for Np in "ACGT-":
                        p_j += prob_t_N(geno, Np) * error_rates[Np][j]
                    prod_j *= p_j
                m_small[g][b] = prod_i * prod_j
        m_big[i] = (m_small * prior_matrix).sum(axis=0)

    p_base_matrix = m_big / (m_big.sum(axis=1))[:, np.newaxis]
    return p_base_matrix


def update_prior_probs():
    global prior_matrix
    likelihood_matrix = np.zeros((genotype_num, H))
    for h in range(H):
        for g, geno in enumerate(all_geno):
            likelihood = likelihood_genotype(geno, alignment[:, h])
            likelihood_matrix[g][h] = likelihood

    geno_data_matrix = likelihood_matrix * prior_matrix
    post_matrix = geno_data_matrix / geno_data_matrix.sum(axis=0)
    prior_matrix = (post_matrix.sum(axis=1) / H)[:, np.newaxis]
    return post_matrix

def update_error_rates():
    global error_rates
    base_prob_matrix = np.zeros((H, K, 5))
    for h in range(H):
        base_prob_matrix[h] = prob_N_data(list(alignment[:, h]))

    base_total_prob = base_prob_matrix.sum(axis=0).sum(axis=0)
    for base in "ACGT":
        for error_base in "ACGT-":
            p_error = 0
            for h in range(H):
                m = base_prob_matrix[h]
                p_error += (m[alignment[:, h]==error_base])[:, bases_index[base]].sum()
            error_rates[base][error_base] = p_error / base_total_prob[bases_index[base]]


# initializition
K, H = alignment.shape
all_geno = all_genotype(2)
genotype_num = len(all_geno)

prior_matrix = np.full(shape=(genotype_num, 1), fill_value=1/genotype_num)
for round_ in range(4):
    tmp = update_prior_probs()
    update_error_rates()
p_m = update_prior_probs()

