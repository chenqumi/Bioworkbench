"""
chenqumi@20181023
"""
import os
import sys

from scipy import stats


def QV2error_rate(Qscore):
    return 10 ** (Qscore/(-10))


def poisson_p(alt, error_num_expect):
    p = 1 - stats.poisson.cdf(alt-1, error_num_expect)
    return p


def vaf_limit(depth, error_rate, significance):
    p = 0.1
    alt = 1
    error_num_expect = depth * error_rate
    while p > significance:
        alt += 1
        p = poisson_p(alt, error_num_expect)
    return alt, alt / depth * 100


def vaf_tendency():
    Qscore = 25
    error_rate = QV2error_rate(Qscore)
    for depth in range(100, 10001, 100):
        alt, vaf = vaf_limit(depth, error_rate, significance)
        print(f"{depth}\t{vaf}")


def main():
    # qiagen lung FFPE panel about 0.15 Mb
    depth = TROUGHPUT / panel_size / sample_num
    Qscore = 25
    error_rate = QV2error_rate(Qscore)
    
    alt, vaf = vaf_limit(depth, error_rate, significance)

    
    with open("lower_detection_limit.log", "w") as O:
        O.write(
            f"sequencing throughput:\t2G\n"
            f"panel size:\t{panel_size/10**6} Mb\n"
            f"sample num:\t{sample_num}\n"
            f"significance:\t{significance}\n"
            f"alt support read num:\t{alt}\n"
            f"coverage:\t{depth}\n"
            f"vaf threshold:\t{vaf}\n"
        )


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print(
            f"\nUsage:{sys.argv[0]} <panel_size[Mb]> <sample_num> <significance>"
        )
        sys.exit()

    panel_size, sample_num, significance = sys.argv[1:4]
    panel_size = float(panel_size) * 10**6
    sample_num = int(sample_num)
    significance = float(significance)


    TROUGHPUT = 2 * 10**9


    main()
    vaf_tendency()
