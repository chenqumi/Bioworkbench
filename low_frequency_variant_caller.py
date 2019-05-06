import numpy as np
from collections import Counter
from scipy.stats import chi2


initial_error_rates = {
    "A": {"A": 0.99,   "C": 0.0025, "G": 0.0025, "T": 0.0025, "-": 0.0025},
    "C": {"A": 0.0025, "C": 0.99,   "G": 0.0025, "T": 0.0025, "-": 0.0025},
    "G": {"A": 0.0025, "C": 0.0025, "G": 0.99,   "T": 0.0025, "-": 0.0025},
    "T": {"A": 0.0025, "C": 0.0025, "G": 0.0025, "T": 0.99,   "-": 0.0025},
    "-": {"A": 0,      "C": 0,      "G": 0,      "T": 0,      "-": 1},
}

bases_index = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "-": 4
}

alignment = np.array([
    ["A", "T", "A", "A", "A", "G", "T"],
    ["A", "T", "A", "A", "A", "G", "T"],
    ["A", "T", "A", "A", "A", "G", "T"],
    ["A", "T", "A", "A", "A", "G", "T"],
    ["A", "T", "A", "A", "A", "G", "T"],
    ["A", "T", "A", "A", "A", "G", "T"],
    ["A", "G", "A", "A", "A", "G", "T"],
    ["A", "G", "A", "A", "A", "G", "T"],
    ["A", "G", "A", "A", "A", "G", "T"],
    ["A", "G", "A", "A", "G", "G", "T"],
    ["A", "G", "C", "A", "C", "G", "T"],
    ["A", "C", "C", "A", "C", "G", "T"],
])


def model_value(likelihood, df, significance):
    """
    v_m = 2*log(L(H_m)) - c_dfm(p)
    """
    v_m = 2 * np.log(likelihood)
    if df > 0:
        v_m -= chi2.cdf(1-significance, df)
    return v_m


def model_frequency(ref, alt, bases_all_reads, f_vector, error_rates):
    """
    str ref
    str alt
    iterable-obj bases_all_reads

    return frequency of alt in bases_all_reads
    """
    # model Mx: NO alt allele
    if ref == alt:
        return np.array([error_rates[ref][eb] for eb in "ACGT-"])
    # model Mxy: two alleles
    mat = np.zeros((K, 5))
    for b, base in enumerate("ACGT-"):
        for i, observed_base in enumerate(bases_all_reads):
            mat[i][b] = error_rates[base][observed_base]

    m2 = (mat * f_vector)[:, [bases_index[alt], bases_index[ref]]]
    alt_ref_f_vector = (m2 / m2.sum(axis=1)[:, np.newaxis]).sum(axis=0) / K

    l = [alt_ref_f_vector[0], ]*5
    l[bases_index[ref]] = alt_ref_f_vector[1]
    updated_f_vector = np.array(l)

    return updated_f_vector


def model_likelihood(f_vector, bases_all_reads):
    likelihood = 1
    for base in bases_all_reads:
        likelihood *= f_vector[bases_index[base]]
    return likelihood
    

def frequency_matrix(refs, f_matrix, error_rates):
    updated_f_matrix = np.zeros((H, 5))
    optimized_alts = ""
    probs_variant = []
    for h in range(H):
        reads = alignment[:, h]
        ref = refs[h]
        possible_alts = Counter(reads).keys()
        v_m_max = -99999
        optimized_f_vector = np.repeat(0.5, 5)
        optimized_alt = ref
        sigma_exp_vm = 0
        for alt in possible_alts:
            df = 0 if ref == alt else 1
            f_vec = model_frequency(ref, alt, reads, f_matrix[h], error_rates)
            likelihood = model_likelihood(f_vec, reads)
            v_m = model_value(likelihood, df, 0.05)
            sigma_exp_vm += np.exp(v_m)
            if v_m > v_m_max:
                v_m_max = v_m
                optimized_f_vector = f_vec
                optimized_alt = alt
        updated_f_matrix[h] = optimized_f_vector
        optimized_alts += optimized_alt
        optimized_exp_vm = np.exp(v_m_max)
        probs_variant.append(optimized_exp_vm / sigma_exp_vm)
    return updated_f_matrix, optimized_alts, probs_variant


# def frequency_matrix(refs, alts, f_matrix, error_rates):
#     updated_f_matrix = np.zeros((H, 5))
#     for h in range(H):
#         ref = refs[h]
#         alt = alts[h]
#         updated_f_vector = model_frequency(ref, alt, alignment[:, h], 
#                                             f_matrix[h], error_rates)
#         updated_f_matrix[h] = updated_f_vector
#     return updated_f_matrix


def update_error_rates(f_matrix, error_rates):
    """
    """
    updated_error_rates = {
        "-": {"A": 0, "C": 0, "G": 0, "T": 0, "-": 1},
    }

    i_error_mat = np.zeros((K, 5))
    all_base_error_mat = np.zeros((H, K, 5))
    for h in range(H):
        for b, base in enumerate("ACGT-"):
            for i, observed_base in enumerate(alignment[:, h]):
                i_error_mat[i][b] = error_rates[base][observed_base]

        m = i_error_mat * f_matrix[h]
        m2 = m / m.sum(axis=1)[:, np.newaxis]
        all_base_error_mat[h] = m2

    base_total_prob = all_base_error_mat.sum(axis=1).sum(axis=0)

    for base in "ACGT":
        for error_base in "ACGT-":
            # p_error: P(base|error_base)
            p_error = 0
            for h in range(H):
                m2_ = all_base_error_mat[h]
                p_error += (
                (m2_[alignment[:, h]==error_base])[:, bases_index[base]].sum())
            updated_error_rates.setdefault(base, {})[error_base] = (
                                p_error / base_total_prob[bases_index[base]])

    return updated_error_rates


def prob_variant(ref, alt, bases_all_reads):
    possible_alts = Counter(bases_all_reads).keys()
    for pa in possible_alts:



K, H = alignment.shape
initial_f_vector = np.repeat(0.5, 5)
initial_f_matrix = np.full(shape=(H, 5), fill_value=0.5)
ROUNDS = 4

updated_f_matrix, alts = frequency_matrix("ATAAAGT", initial_f_matrix,
                                                            initial_error_rates)
updated_error_rates = update_error_rates(updated_f_matrix, initial_error_rates)

for r in range(ROUNDS - 1):
    updated_f_matrix, alts = frequency_matrix("ATAAAGT", updated_f_matrix,
                                                            updated_error_rates)
    updated_error_rates = update_error_rates(updated_f_matrix,
                                                            updated_error_rates)
