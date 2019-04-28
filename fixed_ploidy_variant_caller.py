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


def posterior_prob(genotypes, prior_probs, likelihoods):
    post_probs = {}
    geno_probs = [prior_probs[geno]*likelihoods[geno] for geno in genotypes]
    total_prob = sum(geno_probs)
    i = 0
    for geno in genotypes:
        post_probs[geno] = geno_probs[i] / total_prob
        i += 1

    return post_probs


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
    iterable-obj bases_all_reads, list or str. CANNOT be ndarray
    return a K*5 matrix 

    P(rih = N | data)
    rih means read i at h site
    K is reads num, 5 means 5 possible bases(ACGT-)
    """
    mat_p_bases =np.zeros((len(bases_all_reads), 5))
    for i, actual_base in enumerate(bases_all_reads):
        rest_bases = bases_all_reads[:i] + bases_all_reads[i+1:]

        for b, possible_base in enumerate("ACGT-"):
            
            p_rih_data = 0
            for geno in all_geno:
                tmp = (prior_probs[geno] * 
                       prob_t_N(geno, possible_base) * 
                       error_rates[possible_base][actual_base])
                prod = 1
                for j in rest_bases:
                    p_j = 0
                    for Np in "ACGT-":
                        p_j += prob_t_N(geno, Np) * error_rates[Np][j]
                    prod *= p_j

                p_rih_data += tmp * prod

            p_data = 0
            for geno in all_geno:
                fs = prior_probs[geno]
                prod = 1
                for j in bases_all_reads:
                    p_j = 0
                    for Np in "ACGT-":
                        p_j += prob_t_N(geno, Np) * error_rates[Np][j]
                    prod *= p_j
                p_data += fs * prod

            mat_p_bases[i][b] = p_rih_data / p_data

    return mat_p_bases


def update_fs_er(all_geno, round_num=4):
    """
    iterable-obj all_geno, list or set
    int round_num, [4]
    
    EM for updating fs and error rates.
    fs: prior probs for each genotype
    """
    global prior_probs, error_rates

    for round_ in range(round_num):
        p_geno_data_dict = {} 
        for h in range(H):
            # x = 0
            for g, geno in enumerate(all_geno):
                likelihood = likelihood_genotype(geno, alignment[:, h])
                #likeli_matrix[h][g] = likelihood
                # x += likelihood * prior_probs[geno]
                # P(geno, data) = P(data|geno) * P(geno)
                p_geno_data = likelihood * prior_probs[geno]
                p_geno_data_dict.setdefault(h, {})[geno] = p_geno_data

        sum_p_post = {}
        for h, p_geno_data in p_geno_data_dict.items():
            # P(data) = sum(P(data, genotypes)) = sum(P(data|geno)*P(geno))
            p_data = sum(p_geno_data.values())
            for geno in all_geno:
                # P(geno|data) = P(geno, data) / P(data)
                p_post_geno = p_geno_data[geno] / p_data
                sum_p_post[geno] = sum_p_post.get(geno, 0) + p_post_geno

        # updating prior probs of genotype
        for genotype, v in sum_p_post.items():
            prior_probs[genotype] = v / H

        # updating error rates
        base_prob_matrix = np.zeros((H, K, 5))
        for h in range(H):
            base_prob_matrix[h] = prob_N_data(list(alignment[:, h]))
        base_total_prob = base_prob_matrix.sum(axis=0).sum(axis=0)

        for i, base in enumerate("ACGT"):
            P_base_total = base_total_prob[bases_index[base]]
            # error_bases = "ACGT-"[:i] + "ACGT-"[i+1:]
            for eb in "ACGT-":
                idx = np.argwhere(alignment==eb)
                idx2 = [[x[1], x[0], bases_index[base]] for x in idx]
                P_error = 0
                for j in idx2:
                    P_error += base_prob_matrix[j[0], j[1], j[2]]
                p = P_error / P_base_total
                error_rates[base][eb] = p


# ==============================================================================
# initializition
K, H = alignment.shape
all_geno = all_genotype(2)
initial_prior = 1 / len(all_geno) 
prior_probs = {}
for geno in all_geno:
    prior_probs[geno] = initial_prior

# EM for updating fs and error rates
update_fs_er(all_geno, round_num=4)


# calculate post prob for each genotype
p_geno_data_dict = {} 
for h in range(H):
    # x = 0
    for g, geno in enumerate(all_geno):
        likelihood = likelihood_genotype(geno, alignment[:, h])
        #likeli_matrix[h][g] = likelihood
        # x += likelihood * prior_probs[geno]
        # P(geno, data) = P(data|geno) * P(geno)
        p_geno_data = likelihood * prior_probs[geno]
        p_geno_data_dict.setdefault(h, {})[geno] = p_geno_data

sum_p_post = {}
for h, p_geno_data in p_geno_data_dict.items():
    # P(data) = sum(P(data, genotypes)) = sum(P(data|geno)*P(geno))
    p_data = sum(p_geno_data.values())
    for geno in all_geno:
        # P(geno|data) = P(geno, data) / P(data)
        p_post_geno = p_geno_data[geno] / p_data
        print(geno , p_post_geno)