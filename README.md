# Fixed Ploidy Variant Caller

## 更新时间

2019.5.5

[文档描述](http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/current/index.php?manual=Fixed_Ploidy_Low_Frequency_Detection_tools_detailed_descript.html)

## Caller

Fixed Ploidy Variant Caller使用了Bayesian model和Maximum Likelihood approach算法。通过贝叶斯模型计算出的后验概率来判断一个位点是否是突变。对于任一位点，计算每一种基因型的后验概率，如果所有非纯合reference allele基因型的后验概率之和大于自定义的阈值，那么该位点被认为是突变位点，后验概率最大的基因型即为该位点的基因型。

用最大似然估计(maximum likelihood estimate)来估计模型中的未知参数，由于该问题中有两个未知参数，所以使用最大期望(Expectation Maximization)方法来求解。

## Model

模型由**the possible site types, S**, **their prior probabilities, $f_s, s \in S​$**和**the sequencing errors, e**组成

### Prior site type probabilities

可能的基因型完全由倍性(ploidy)决定，对于二倍体，$$S=\{AA, AC, AG, AT, A-, CC, CG, CT, C-, GG, GT, G-, TT, T-, --\}$$

将$f_s$记为每一基因型的先验概率，是比对数据中真实位点类型的频率(frequencies of the true site types in the mapping)，是一个未知的值。

### Error probabilities

sequencing errors模型描述了测序仪将真实碱基N测为碱基M的概率，$ \{e(N \to M)|N,M \in \{A,C,G,T,-\}\}$，该值也是未知的。

### Posterior probabilities

$P(t|data)\ =\ \frac{P(data|t)P(t)} {P(data)}​$

$= \frac{P(data|t)P(t)} {\sum_\limits{s \in S} P(data|s)P(s)}$ (式1.)

其中：

$P(t)​$：基因型为t的先验概率

data：覆盖位点的所有碱基。假设该位点有k条read覆盖，$P(data|t) = P(n_1, ..., n_k | t)​$



记$P_t(N)$为基因型为t时其中一个等位基因为N的概率。$P_t(N)$的值也是由倍性决定的，如：$P_{AA}(A) = 1$, $P_{AC}(A) = 0.5$, $P_{AC}(G) = 0$, $P_{ACC}(A) = \frac{1}{3}$

基于上述定义，我们可以将$P(n_1, ..., n_k | t)$写为：

$P(n_1, ..., n_k | t) \ =\ P(n_1 | t)...P(n_k | t)$

$= \prod_\limits{i=1}^{k} P(n_i | t)​$

$= \prod_\limits{i=1}^{k} \sum_\limits{N \in \{A,C,G,T,-\}} P_t(N) \times e_q(N \rightarrow n_i)​$  (式2.)



将式2.代入式1.

$P(t | n_1, ..., n_k) = \frac{P(n_1, ..., n_k | t) f_t} {\sum_\limits{s \in S} P(n_1, ..., n_k | s) f_s}​$

$= \frac{f_t \prod_\limits{i=1}^{k} \sum_\limits{N \in \{A,C,G,T,-\}} P_t(N) \times e_q(N \rightarrow n_i)} {\sum_\limits{s \in S} f_s \prod_\limits{i=1}^{k} \sum_\limits{N \in \{A,C,G,T,-\}} P_s(N) \times e_q(N \rightarrow n_i)}​$ (式3.)

## Estimating parameters

### Updating prior probabilities

在式3.中，我们得到的是某一位点的后验概率，现考虑alignment中的所有位点。"Let $h$ index the sites in the alignment (h = 1, .., H)"

$f_{t}^{*} = \frac{\sum_\limits{h=1}^{H} \frac{f_t \prod_\limits{i=1}^{k} \sum_\limits{N \in \{A,C,G,T,-\}} P_t(N) \times e_q(N \rightarrow n_i)} {\sum_\limits{s \in S} f_s \prod_\limits{i=1}^{k} \sum_\limits{N \in \{A,C,G,T,-\}} P_s(N) \times e_q(N \rightarrow n_i)} } {H}​$ (式4.)

### Updating error probabilities

对于一个特定位置$r_{i}^{h}$, read i, site h, 它在测序数据为$n_{1}^{h}, ..., n_{k}^{h}$的情况下，真实碱基为N的联合概率为：

$P({r_i^h}=N, {n_1^h}, ..., {n_k^h})$ 联合概率

$= \sum_\limits{s \in S} f_s P({r_i^h}=N, {n_1^h}, ..., {n_k^h} | s)$ 写成全概率的形式

$= \sum_\limits{s \in S} f_s P({r_i^h}=N, {n_i^h} | s) P({n_1^h | s})...P({n_{i-1}^h | s})P({n_{i+1}^h | s})...P({n_k^h | s})$ 拆分全概率中联合概率为独立事件概率的乘积

$= \sum_\limits{s \in S} f_s P({r_i^h}=N, {n_i^h} | s) \prod_\limits{j \ne i} P({n_j^h} | s)$ 简写上一式

$= \sum_\limits{s \in S} f_s (P_s(N) \times e_{q_i^h}(N \rightarrow {n_i^h}) \prod_\limits{j \ne i} \sum_\limits{N' \in \{A,C,G,T,-\}} P_s(N') \times e_{q_j}({N'} \rightarrow {n_j^h}))$  (式5.)



运用贝叶斯公式：

$P({r_i^h}=N | {n_1^h}, ..., {n_k^h}) = \frac{P({r_i^h}=N, {n_1^h}, ..., {n_k^h})} {P({n_1^h}, ..., {n_k^h})}$

$ = \frac{P({r_i^h}=N, {n_1^h}, ..., {n_k^h})} {\sum_\limits{N' \in \{A,C,G,T,-\}} P({r_i^h}={N'}, {n_1^h}, ..., {n_k^h})}$ (式6.)



 将 式5. 代入 式6.

$P({r_i^h}=N | {n_1^h}, ..., {n_k^h})$

$= \frac{\sum_\limits{s \in S} f_s (P_s(N) \times e_{q_i^h}(N \rightarrow {n_i^h}) \prod_\limits{j \ne i} \sum_\limits{N' \in \{A,C,G,T,-\}} P_s(N') \times e_{q_j}({N'} \rightarrow {n_j^h}))} {\sum\limits_{{N' \in \{A,C,G,T,-\}}} (\sum_\limits{s \in S} f_s (P_s(N') \times e_{q_i^h}(N' \rightarrow {n_i^h}) \prod_\limits{j \ne i} \sum_\limits{N'' \in \{A,C,G,T,-\}} P_s(N'') \times e_{q_j}({N''} \rightarrow {n_j^h})))}$

$= \frac{\sum_\limits{s \in S} f_s (P_s(N) \times e_{q_i^h}(N \rightarrow {n_i^h}) \prod_\limits{j \ne i} \sum_\limits{N' \in \{A,C,G,T,-\}} P_s(N') \times e_{q_j}({N'} \rightarrow {n_j^h}))} {\sum_\limits{s \in S} f_s \prod_\limits{j} \sum_\limits{N'' \in \{A,C,G,T,-\}} P_s(N'') \times e_{q_j^h}(N'' \rightarrow {n_j^h})}$ (式7.)

式7.计算了特定位置碱基为N的概率，因此我们可以计算一次alignment中所有位置的概率以更新错误率：

$e_{q}^{*}(N \rightarrow M) = \frac{\sum_\limits{h} \sum_\limits{i=1, .., k: n_{i}^{h}=M} P(r_{i}^{k} = N | n_{1}^{k}, ..., n_{k}^{h})} {\sum_\limits{h} \sum_\limits{i=1, .., k} P(r_{i}^{k} = N | n_{1}^{k}, ..., n_{k}^{h})}$





