# MOLECULAR ECOLOGY RESOURCE | 软件介绍--PIXY ：对核苷酸多样性(π)和DXY的无偏估计

![文献来源](https://gitee.com/masterttz/we-chat/raw/master/img/20220215161716.png)

![软件说明文档：https://pixy.readthedocs.io/en/latest/](https://gitee.com/masterttz/we-chat/raw/master/img/20220215161834.png)

## 什么是PIXY？

pixy 是一个命令行工具，用于从 VCF 文件计算群体的核苷酸多样性pi**（每个位点的平均杂合度）和**dxy（每个位点的群体之间的平均核苷酸差异数）。

许多用于从 VCF 计算 pi 和 dxy 的工具在存在缺失数据的情况下会产生有偏差的估计。这是因为这些方法通常会做出简化假设，即如果基因型缺失，则按状态是纯合参考 (0/0)。有关详细信息，请参阅[pixy 的论文](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13326)

![Figure 1 from Korunes & Samuk 2021](https://gitee.com/masterttz/we-chat/raw/master/img/men13326-fig-0001-m.jpg)

> In Case 1, all missing data is assumed to be present but invariant. This results in a deflated estimate of π. In Case 2, missing data are simply omitted from the calculation, both in terms of the number of sites (the final denominator) and the component denominators for each site (the n choose 2 terms). This results in an unbiased estimate of π. The adjusted π method (Case 2) is implemented for VCFs in pixy. Invariant sites are represented as sites with no ALT allele, and greyed-out sites are those that failed to pass a genotype filter requiring a minimum number of reads covering the site (Depth>=10 in this case).

## pixy的优势

- 快速高效地处理非变异位点 VCF
- 计算任意数量种群的π 和 dxy
- 在任意大小的窗口中计算所有统计数据，输出包含所有计算的所有原始信息。
- 这大部分是通过广泛使用的 python 库 scikit-allel 中的现有数据结构和函数来实现的

# 安装

## 通过conda

pixy 目前可通过 conda-forge 安装在 Linux/OSX 系统上。您可以使用以下两个 conda install 命令安装它：

```
conda install -c conda-forge pixy
conda install -c bioconda htslib
```

您可以通过下列来测试安装是否成功：

```
pixy --help
```

- 有关安装 conda 的信息：

  [anaconda（更多功能和初始模块）](https://docs.anaconda.com/anaconda/install/) 

  [miniconda（安装更快）](https://docs.conda.io/en/latest/miniconda.html)

## 通过 Git 或下载

或者，您可以从 Github 克隆，[手动安装依赖项](https://github.com/ksamuk/pixy)

# 参数命令

下面是 pixy 接受的必需和可选参数的列表。

## 必需的：

- **--stats [pi|fst|dxy]**

  从 VCF 计算哪些统计数据（pi、dxy 和/或 fst，用空格分隔）

- **--vcf [路径/到/vcf.vcf.gz]**

  输入 VCF 的路径（bgzipped 和 tabix 索引）。

- **--populations [路径/到/populations_file.txt]**

  群体分类信息文件的路径。

## 此外还需包含：

- **--window_size [整数]**

  计算 pi/dxy 的碱基对窗口大小。

- **--bed_file [路径/到/populations_file.txt]**

  包含计算汇总统计信息的区域（chrom、chromStart、chromEnd）的BED 文件的路径

## 可选参数：

- **--n_cores[整数]**

  用于并行处理的 CPU 数量（默认值 = 1）。

- **--output_folder [路径/到/输出/文件夹]**

  将写入输出的文件夹，例如 path/to/output_folder，默认为当前工作目录。

- **--output_prefix [前缀]**

  输出文件的可选前缀，例如“输出”将导致写入[输出文件夹]/output_pi.txt，默认为“pixy”。

- **--chromosomes ['list,of,chromosomes']**

  单引号、逗号分隔的染色体列表（例如，'X,1,2'）。默认为 VCF 中的所有染色体。

- **--interval_start [整数]**

  计算 pi/dxy 的特定间隔的开始。仅在计算单个染色体时有效。默认为 1。

- **--interval_end [整数]**

  计算 pi/dxy 的特定间隔的结束。仅在计算单个染色体时有效。默认为染色体的最后一个位置。

- **--sites_file [路径/到/sites_file.txt]**

  包含要（仅）包含在计算中的无标题站点列表（CHROM、POS）的制表符分隔文件的路径

- **--chunk_size [整数]**

  在任何给定时间从 VCF 读取的大约站点数。默认为 100000。较小的数字可以减少内存使用。

- **--fst_type [wc|哈德森]**

  要使用的 FST 估计器，其中之一是：“wc”（Weir and Cockerham 1984）或“hudson”（Hudson 1992, Bhatia et al. 2013）。默认为“wc”

- **--bypass_invariant_check [否|是]**

  忽略对非变异位点的检查。谨慎使用！

- **- version**

  打印软件版本号。

- **--citation**

  打印 pixy 的参考文献。

- **- help**

  打印帮助信息。

- **- silent**

  禁止所有输出。

### 一个例子：

```css
pixy --stats pi fst dxy \
--vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--window_size 10000 \
--n_cores 4 \
--output_folder output \
--output_prefix pixy_output
```

# 输入文件

## 群体分类

*描述*

人口文件用于将个人分配给人口。对于成对统计（dxy 和 fst），将生成总体之间的所有成对比较。它在命令行中使用`--populations`参数指定。

*格式*

无标题、制表符分隔的文件，第一列包含样本名称（如 VCF 中所示），第二列包含任意群体 ID。

*例子*

```
ERS223827   BFS
ERS223759   BFS
ERS223750   BFS
ERS223967   AFS
ERS223970   AFS
ERS223924   AFS
ERS224300   AFS
ERS224168   KES
ERS224314   KES
```

## BED

*描述*

BED 文件用于手动定义用于计算汇总统计的基因组区域。它在命令行中以`--bed_file`.

请注意，pixy 期望间隔是一个索引，*而不是*像许多典型的 BED 文件那样索引为零。如果您的 BED 文件是零索引的，您需要转换间隔（例如，通过在开始和结束处添加 1）。

*格式*

无标题、制表符分隔的文件，第一列包含染色体 ID，第二列包含窗口的第一个位置，第三列包含窗口的最后一部分。每行 = 一个窗口。

*例子*

```css
X   1       10000
X   10001   20000
X   20001   30000
X   30001   40000
X   40001   50000
X   50001   60000
1   1       10000
1   10001   20000
1   20001   30000
1   30001   40000
1   40001   50000
1   50001   60000
```

## 位点文件

*描述*

位点用于将计算限制在特定位点（例如 4 倍简并位点、内含子等）。它在命令行中以`--sites_file`.

*格式*

无标题、制表符分隔的文件，第一列包含染色体 ID (CHROM)，第二列包含单个位点的位置 (POS)。每行 = 一个站点。

*例子*

```css
X   1
X   2
X   3
X   4
X   5
X   6
X   7
X   8
X   9
X   10
X   11
X   12
X   13
X   14
X   15
1   1
1   2
1   3
1   4
1   5
1   6
1   7
1   8
1   9
1   10
1   11
1   12
1   13
1   14
1   15
```

# 生成包含非变异位点的VCF

pixy 通过允许用户以标准文件格式输入数据来促进 π 和 dxy 的正确计算：变异文件格式 (VCF)。AllSites VCF 包含除了变异位点之外的非变异位点。一些常用的软件可以生成 AllSites VCF。

在这里，我们提供了使用两种最广泛使用的工具生成 AllSites VCF 的快速入门指南：BCFtools 和 GATK。只要给定一组对齐的数据（BAM 文件）和用于对齐它们的参考序列，就可以运行这些工具中的任何一个。

**BCFtools 和 GATK 也可以很好地过滤 VCF，我们建议在使用 pixy 分析之前利用这一点过滤您的数据。**

> **利用基因组间隔改进运行时间：**如果生成 AllSites VCF 非常耗时，我们建议通过将分析分解为更小的基因组区域来并行化您的流程。在我们的测试数据集中，我们并行运行单个染色体。根据数据集中染色体的大小，将染色体分解成更小的间隔以进行call SNP可能是有益的。可以使用 GATK 中的 -L 参数或 bcftools mpileup 中的 -r 参数指定基因组间隔。

## 使用 BCFtools (mpileup/call) 生成 AllSites VCF

BCFtools mpileup 可用于产生基因型似然性，此操作后可调用 bcftools 调用 SNPs/INDELS。BCFtools 提供了许多灵活的选项，记录在这里：[https ://samtools.github.io/bcftools/bcftools-man.html#call](https://samtools.github.io/bcftools/bcftools-man.html#call)

在此示例中，我们调用 mpileup 并将输出传递给调用变体并生成 AllSites VCF：

```
bcftools mpileup -f <reference.fa> -b <bamlist.txt> -r <X> | bcftools call -m -Oz -f GQ -o <output>
```

此处选择的选项的注意事项：

- b 将 mpileup 指向文本文件中包含的 BAM 文件列表。
- r 指定基因组区域。在本例中，我们指定 X 染色体，但 mpileup 提供了多种用于指定区域的选项（CHR|CHR:POS|CHR:FROM-TO|CHR:FROM-[,…]）。或者， -R <file> 可用于从提供的文件中读取区域。
- Oz 指定压缩的 VCF 输出
- f GQ 表示应该为每个样本输出 FORMAT 字段 GQ

## 使用 GATK 生成 AllSites VCF

GATK 建议首先在 GVCF 模式下使用 HaplotypeCaller 调用每个样本的变体（下面的步骤 1）。接下来，GenomicsDBImport 整合来自 GVCF 文件的跨样本信息，以提高联合基因分型的效率（下面的第 2 步）。在第 3 步中，GenotypeGVCFs 生成一组联合称为 SNP 和 INDELS 的集合，可用于过滤和分析。我们建议您查阅此处的完整 GATK 文档：[https ://software.broadinstitute.org/gatk/](https://software.broadinstitute.org/gatk/)

第 1 步 - HaplotypeCaller（应为每个 BAM 文件运行此步骤）：

```css
gatk-4.0.7.0/gatk --java-options "-Xmx4G" HaplotypeCaller \
-R <reference.fa> -I <input.bam> -O <output.g.vcf> -ERC GVCF -L <X>
```

第 2 步 - 基因组学DBImport：

```css
gatk-4.0.7.0/gatk --java-options "-Xmx4g" GenomicsDBImport \
-V $file1 -V $file2 --genomicsdb-workspace-path <allsamples_genomicsdb> \
-L <X>
```

第 3 步 - 基因型 GVCF：

```css
gatk-4.0.7.0/gatk --java-options "-Xmx4g" GenotypeGVCFs \
-R <reference.fa> -V gendb://<allsamples_genomicsdb> \
-all-sites -L <X> -O <output-allsites.vcf.gz>
```

关于上面选择的选项的注释：

- ERC GVCF 设置发出参考置信度分数的模式。GVCF 格式指定压缩的非变体块。
- L 指定基因组区域。在本例中，我们指定 X 染色体
- V 指定输入的变量数据。对于 GenomicsDBImport，这是 GVCF 文件。对于 GenotypeGVCFs，联合基因分型的变异数据作为 GenomicsDB 工作区提供，该工作区是在上一步中使用 GenomicsDBImport 创建的。
- all-sites 表明最终的 VCF 应包括基因分型后发现为非变异的基因座。对我们来说最重要的参数。

重要提示：在 GATK v4.x 中，我们发现在运行 GenotypeGVCF 时必须指定间隔 (-L)。如果没有指定的间隔，它似乎会遇到缺少参考置信块，导致它失败。即使有问题的块在传递给它的 GenomicsDB 区间之外也是如此。我们建议始终指定一个间隔 (-L) 以避免此类问题。

> **GATK 版本：**在上面的示例中，我们使用 GATK v4，但在 GATK v3.x 中也可以通过运行带有“-allSites”参数的 GenotypeGVCFs 轻松生成 AllSites VCF。（请注意与 GATK v4 中的“-all-sites”的语法略有不同）。

# 分步安装和使用

笔记

pixy 目前仅适用于 Linux 和 macOS 系统。

## 1. 生成具有 Invariant Sites 的 VCF 并执行过滤

如果您尚未生成“allsites”VCF（具有不变站点的 VCF），[请参阅此处的指南](https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html)。

笔记

在处理全基因组数据时，我们建议您*为每个染色体生成单独的不变位点 VCF*。这是为了防止内存限制问题。自然，这对于减少表示排序来说不是什么问题。

我们建议使用专用于在 VCF 上执行此类操作的标准工具：VCFtools 和 BCFtools（均在 bioconda 上提供）。

### 站点级过滤

对种群遗传学的 VCF 过滤进行全面处理超出了我们的范围，但作为一个好的开始，请参阅：[https](https://speciationgenomics.github.io/filtering_vcfs/) ://speciationgenomics.github.io/filtering_vcfs/ 。

在最低限度，我们建议 (a) 使用 GATK 最佳实践站点过滤程序（如果适用）过滤您的 VCF，并且 (b) 对缺失、站点质量得分和平均深度（最小和最大）执行进一步过滤。

这是一个使用 VCFtools 的示例。具体值（尤其是 min/max-meanDP）将根据您的数据集而有所不同：

```css
vcftools --gzvcf my_vcf.vcf.gz \
--remove-indels \
--max-missing 0.8 \
--min-meanDP 20 \
--max-meanDP 500 \
--recode --stdout | gzip -c > my_filtered_vcf.vcf.gz
```

> 从 GATK 版本 4.2.0.0 开始，GenotypeGVCFs 仅在没有缺失基因型的情况下将 QUAL 分数分配给不变位点。如果您希望在 QUAL (-minQ) 上进行过滤，则需要分别过滤不变量和变体站点**，**而 QUAL 过滤器仅适用于变体站点（有关详细信息，请参见下文）。

### 可选：种群遗传过滤器

根据您的目标，您还可以考虑过滤掉具有严重 HWE 违规的站点（尝试使用 VCFtools 使用 –hwe 0.001）、观察到的异常高的杂合性或等位基因深度不平衡。有关这些注意事项的更多详细信息，请参阅本文https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12613。

如果您的研究生物具有高水平的基因复制（例如，在许多植物和鱼类物种中，在全基因组复制后重新二倍化），最后两个考虑因素尤其重要。

如果您的 VCF 包含变异位点和不变位点（此时应该如此），应用基于群体遗传的过滤器将导致您的不变位点丢失。为避免这种情况，请分别过滤不变和变体站点并将两个结果文件连接起来。以下是使用 VCFtool 和 BCFtools 实现此目的的一种方法示例：

> ```css
> #!/bin/bash
> # requires bcftools/bgzip/tabix and vcftools
> 
> # create a filtered VCF containing only invariant sites
> vcftools --gzvcf test.vcf.gz \
> --max-maf 0 \
> [add other filters for invariant sites here] \
> --recode --stdout | bgzip -c > test_invariant.vcf.gz
> 
> # create a filtered VCF containing only variant sites
> vcftools --gzvcf test.vcf.gz \
> --mac 1 \
> [add other filters for variant sites here] \
> --recode --stdout | bgzip -c > test_variant.vcf.gz
> 
> # index both vcfs using tabix
> tabix test_invariant.vcf.gz
> tabix test_variant.vcf.gz
> 
> # combine the two VCFs using bcftools concat
> bcftools concat \
> --allow-overlaps \
> test_variant.vcf.gz test_invariant.vcf.gz \
> -O z -o test_filtered.vcf.gz
> ```

## 2.安装conda

如果您还没有安装 Anaconda https://docs.anaconda.com/anaconda/install/

## 3. 创造新环境

为使用 pixy 创建并激活一个新的 conda 环境：

```css
conda create --name pixy
conda activate pixy
```

## 4.安装PIXY

通过 conda-forge 频道安装 pixy。还要从 bioconda 安装所需的 htslib 包。

```css
conda install --yes -c conda-forge pixy
conda install --yes -c bioconda htslib
```

要查看参数列表并测试 pixy 安装，请键入：

```css
pixy --help
```

## 5. 创建群体文件

创建人口文件。这是一个无标题、制表符分隔的文件，其中第一列包含样本名称（与 VCF 中所表示的完全相同），第二列包含总体名称（可以是任何名称，但应该保持一致！）。

例如：

```css
ERS223827   BFS
ERS223759   BFS
ERS223750   BFS
ERS223967   AFS
ERS223970   AFS
ERS223924   AFS
ERS224300   AFS
ERS224168   KES
ERS224314   KES
```

## 6.运行小精灵

运行pixy！一个例子如下所示。

```css
pixy --stats pi fst dxy \
--vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
--populations data/vcf/ag1000/Ag1000_sampleIDs_popfile.txt \
--window_size 10000 \
--n_cores 4 \
--chromosomes 'X'
```

> pixy 会忽略非双等位基因位点和 INDEL，即使它们在预过滤后仍留在 VCF 中。

## 7. 利润

解析输出文件，享受您对 pi 和 dxy 的无偏估计！

## 8. 保持最新

您可以通过重新运行使 pixy 保持最新状态：

```css
conda install --yes -c conda-forge pixy
```

您可以通过以下方式检查您是否拥有最新版本：

> ```css
> pixy --version
> ```

# 使用示例

## 基本用法

```
pixy --stats pi fst dxy \
--vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--window_size 10000 \
--n_cores 2
```

## 专注于特定的基因组区间

```
pixy --stats pi fst dxy \
--vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--chromosomes 'X' \
--interval_start 15000 \
--interval_end 400000 \
--window_size 10000
```

## 使用 .BED 文件定义窗口

当使用 BED 文件定义时，窗口可以是任意大小且不均匀。

```
pixy --stats pi fst dxy \
--vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--bed_file genomic_windows.bed
```

## 使用站点文件将计算限制到特定站点

站点文件定义将包含在计算中的各个站点。这可用于获得 4 倍简并 pi，或特定类别的基因组元素（基因、内含子等）的统计数据。

```
pixy --stats pi fst dxy \
--vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--sites_file valid_sites.txt
```

## 站点级估算

注意：站点级估计的计算速度将比窗口估计慢得多。

```
pixy --stats pi fst dxy \
--vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--window_size 1 \
--chromosomes 'X' \
--interval_start 15000 \
--interval_end 400000
```

[下一个](https://pixy.readthedocs.io/en/latest/output.html)[以前的](https://pixy.readthedocs.io/en/latest/guide/pixy_guide.html)



## 参考文献

Korunes, K.L. and Samuk, K. (2021), pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Molecular Ecology Resources. Accepted Author Manuscript. https://doi.org/10.1111/1755-0998.13326

