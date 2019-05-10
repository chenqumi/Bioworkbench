## 关于QCI-A中coverage的cut-off值设定问题

### 更新时间

2018.10.23

### 问题描述

在分析时发现一个位点的覆盖度可以达到1000以上，但过滤的阈值仅为200，设置的依据是什么

### 理解

1. GeneReader的通量为2G，Lung FFPE panel的大小约为0.15M，一次测序上4个样，则平均覆盖度

   avg_coverage = throughput / panel_size / sample_num 

   约为3333X

   区域的覆盖度并不均一，较低的阈值可以保留一些区域。阈值设置过高会遗漏那些低覆盖的true positive位点

2. 提高测序深度能够下探variant的检测下限，当significance为1e-5时，100深度的位点能够检测的vaf约为4-6%，和QCI-A参数设置基本一致。应该是综合上样量、检测准确度等因素做出的综合考虑