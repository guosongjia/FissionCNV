# Breakpoint annotation columns (res/final/{sample}.bed)

`res/final/{sample}.bed` 在启用 `breakpoint-analysis` 后包含 38 列：原有 16 列 CNV 信息 + 22 列断点注释。

## 原有列（1-16，由 merge → duphold → filter → annotate 产生）

| 列名 | 说明 |
|------|------|
| chromosome | 染色体 |
| start | CNV 起始坐标（0-based） |
| end | CNV 终止坐标 |
| cn | 拷贝数（haploid） |
| cnv | deletion 或 duplication |
| accumScore | 各工具置信分数累计 |
| dupholdScore | duphold 综合质量分（0-100） |
| dhfc, dhbfc, dhffc | duphold 三个 fold-change 指标 |
| toolName | 支持该 CNV 的工具列表（逗号分隔） |
| toolNum | 支持工具数量 |
| gc | CNV 区域 GC 含量 |
| sample | 样本名 |
| GS | low-complexity 评分（100 = 无重叠） |
| MS | mappability 评分（100 = 无重叠） |

## 断点注释列（17-38）

精确工具 CNV（toolName 含 smoove/delly）：所有列都有值。
Depth-only CNV（仅 cnvkit/cnvpytor/mops）：mh 列为空，CI 列为 0，repeat 和 BAM 证据列正常填写。

| 列名 | 说明 |
|------|------|
| ci_5p_left | 5' 断点向左的置信区间（bp），来自 smoove/delly CIPOS（depth-only 为 0） |
| ci_5p_right | 5' 断点向右的置信区间（bp） |
| ci_3p_left | 3' 断点向左的置信区间（bp），来自 smoove/delly CIEND |
| ci_3p_right | 3' 断点向右的置信区间（bp） |
| strict_mh_len | 严格微同源性长度（0 = NHEJ，≥4 = MMEJ 候选）；depth-only 留空 |
| fuzzy_mh_len | 模糊微同源性长度（允许 1 个错配） |
| mh_seq | 微同源序列（strict 部分的碱基） |
| dup_assumption | DUP 微同源性计算假设（tandem）；DEL/depth-only 留空 |
| in_repeat | 任一断点窗口命中已知重复元件 |
| repeat_type_5p | 5' 断点窗口命中的重复元件（按优先级取一个） |
| repeat_type_3p | 3' 断点窗口命中的重复元件 |
| repeat_types_all_5p | 5' 断点窗口命中的全部重复元件类型（逗号分隔） |
| repeat_types_all_3p | 3' 断点窗口命中的全部重复元件类型 |
| simple_overlap_5p | 5' 断点窗口是否落在低复杂度区域（longdust） |
| simple_overlap_3p | 3' 断点窗口是否落在低复杂度区域 |
| split_reads_5p | 5' 断点 ±10bp 内带 SA tag 的 read 数（结构变异证据） |
| split_reads_3p | 3' 断点 ±10bp 内带 SA tag 的 read 数 |
| discordant_pairs_5p | 5' 断点 ±10bp 内异常 read pair 数 |
| discordant_pairs_3p | 3' 断点 ±10bp 内异常 read pair 数 |
| nested | 该 CNV 是否完全嵌套在同样本的另一个 CNV 内 |
| at_chrom_boundary | 断点窗口是否触及染色体边界 |
| dup_too_short | DUP 长度 <50bp 时跳过微同源性计算（depth-only 留空） |

## 设计说明

- **微同源性仅对精确工具 CNV 计算**：smoove/delly 给出的断点坐标在多数情况下精度足够算 mh；depth 工具的 bin 级边界精度（kb 级）不适合 bp 级 mh 分析。
- **BAM 证据对所有 CNV 都计算**：split-read 和 discordant pair 计数可以独立判断 CNV 是否有结构变异支持，即使 caller 没给精确坐标。
- **重复元件优先级**：KMDs > LTR > tRNA > rRNA > centromeric > subtelomeric > none（命名小元件优先于大区域类别）。
- **Repeat 注释文件均可选**：`config.yaml` 的 `data.breakpoint-annotations` 各字段留空时，对应类别会自动跳过。
