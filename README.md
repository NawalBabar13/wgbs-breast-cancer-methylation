#  WGBS — Whole Genome Bisulfite Sequencing Analysis of Breast Cancer Methylomes

##  Assignment Overview
This project performs a complete **Whole Genome Bisulfite Sequencing (WGBS)** analysis pipeline on breast cancer samples using the **Galaxy bioinformatics platform**. The analysis reproduces and extends the methodology from:

> Lin et al. (2015) *"Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes"* — PLOS ONE. DOI: 10.1371/journal.pone.0118453

---

##  Background & Theory

### What is DNA Methylation?
DNA methylation is an **epigenetic modification** where a methyl group (CH₃) is added to cytosine bases at **CpG dinucleotide sites**. It controls gene expression without changing the DNA sequence:
- **Methylated CpG** → gene is **silenced**
- **Unmethylated CpG** → gene is **active**

### Why Does Methylation Matter in Cancer?
In cancer, **tumor suppressor genes** become hypermethylated and silenced, allowing uncontrolled cell growth. Meanwhile other regions become hypomethylated, causing genomic instability. Studying these changes helps identify dysregulated genes in cancer.

### What is WGBS?
WGBS maps DNA methylation at **single-base resolution** across the entire genome using **bisulfite treatment**:
1. Bisulfite converts **unmethylated C → T**
2. **Methylated C stays as C**
3. Comparing C vs T at each CpG site reveals methylation status

This is why QC reports show unusually high T% and low C% — that is **expected and correct**!

### Study Samples (Lin et al. 2015)
| Sample | Type | Description |
|--------|------|-------------|
| NB | Normal Breast | Healthy normal breast tissue |
| BT089 | Benign Tumor | Fibroadenoma |
| BT126 | Malignant Tumor | Invasive ductal carcinoma |
| BT198 | Malignant Tumor | Invasive ductal carcinoma |
| MCF7 | Cancer Cell Line | Breast adenocarcinoma |

---

##  Tools Used

| Tool | Version | Purpose |
|------|---------|---------|
| Galaxy | usegalaxy.org | Cloud bioinformatics platform |
| Falco | 1.2.4+galaxy0 | Quality control of raw reads |
| bwameth | 0.2.7+galaxy0 | Bisulfite-aware alignment |
| MethylDackel | 0.5.2+galaxy0 | Methylation extraction |
| Wig/BedGraph-to-bigWig | 3.5.4+galaxy0 | Format conversion |
| computeMatrix | 3.5.4+galaxy0 | Signal computation around features |
| plotProfile | 3.5.4+galaxy0 | Methylation profile visualization |
| Metilene | 0.2.6+galaxy0 | DMR identification |

**Reference Genome:** Human hg38 (GRCh38)

---

##  Complete Pipeline

---

###  Step 1: Data Upload

Uploaded paired-end bisulfite sequencing FASTQ files into Galaxy:
https://zenodo.org/record/557099/files/subset_1.fastq
https://zenodo.org/record/557099/files/subset_2.fastq

- 458,003 sequences per file
- Read length: 20–150 bp
- Encoding: Sanger / Illumina 1.9

---

### Step 2: Quality Control with Falco

Ran Falco on both raw FASTQ files to assess read quality.

**Tool:** Falco v1.2.4+galaxy0
**Input:** subset_1.fastq + subset_2.fastq
**Output:** HTML quality reports

#### Dataset 1 (subset_1.fastq) Results:

**Basic Statistics:**
![Falco Dataset 1 Basic Stats](falco_dataset1_basic_stats.png)

**Per Base Sequence Quality:**
![Falco Dataset 1 Per Base Quality](falco_dataset1_per_base_quality.png)

**Per Sequence Quality Scores:**
![Falco Dataset 1 Per Seq Quality](falco_dataset1_per_seq_quality.png)

**Per Base Sequence Content:**
![Falco Dataset 1 Sequence Content](falco_dataset1_seq_content.png)

**Per Sequence GC Content:**
![Falco Dataset 1 GC Content](falco_dataset1_gc_content.png)

**Per Base N Content:**
![Falco Dataset 1 N Content](falco_dataset1_n_content.png)

**Sequence Length Distribution:**
![Falco Dataset 1 Seq Length](falco_dataset1_seq_length.png)

**Sequence Duplication Levels:**
![Falco Dataset 1 Duplication](falco_dataset1_seq_duplication.png)

**Adapter Content:**
![Falco Dataset 1 Adapter](falco_dataset1_adapter_content.png)

#### Dataset 2 (subset_2.fastq) Results:

**Full QC Report:**
![Falco Dataset 2](falco_dataset2_basic_stats.png)

| Metric | Result | Note |
|--------|--------|------|
| Basic Statistics | ✅ PASS | 458,003 sequences |
| Per base sequence quality | ✅ PASS | Quality >30 |
| Per sequence quality scores | ✅ PASS | Sharp high-quality peak |
| Per base sequence content | ❌ FAIL | **Expected!** Bisulfite C→T conversion |
| Per sequence GC content | ❌ FAIL | **Expected!** GC skewed by conversion |

> The FAIL results are **NOT errors** — they are the expected signature of successful bisulfite conversion!

---

###  Step 3: Alignment with bwameth

Aligned bisulfite reads to hg38 using bwameth — a bisulfite-aware aligner.

**Why not a normal aligner?** Normal aligners reject bisulfite reads because converted C→T look like mismatches. bwameth handles this correctly.

**Tool:** bwameth v0.2.7+galaxy0
**Input:** subset_1.fastq + subset_2.fastq + hg38 reference
**Output:** BAM file (119.9 MB)

**IGV Visualization (chr1:5,190,700-5,191,200):**
![IGV Alignment](igv_visualization.png)

Individual paired-end reads are visible aligned to chromosome 1. Colors indicate nucleotide mismatches at each position.

---

###  Step 4: Methylation Extraction with MethylDackel

Extracted CpG methylation fractions from the aligned BAM file.

**Tool:** MethylDackel v0.5.2+galaxy0
**Input:** BAM file from bwameth
**Output:** bedGraph file (~607,800 CpG sites, 17.6 MB)

**Output Preview:**
![MethylDackel Output](methyldackel_cpg.png)

**bedGraph format:**
| Column | Meaning |
|--------|---------|
| Chr | Chromosome |
| Start | CpG start position |
| End | CpG end position |
| Value | Methylation fraction (0.0–1.0) |

---

###  Step 5: Methylation Profile Visualization

Used computeMatrix + plotProfile to visualize methylation around CpG islands and TSS.

**Single Sample Profile:**
![plotProfile Single Sample](plotprofile_single_sample.png)

**All Breast Cancer Samples Comparison:**
![plotProfile All Samples](plotprofile_all_samples.png)

The profiles show methylation levels from -1kb to +1kb around CpG islands. Normal breast (NB) shows distinct patterns compared to cancer samples (BT089, BT126, BT198, MCF7).

---

###  Step 6: Differentially Methylated Regions with Metilene

Identified genomic regions significantly differently methylated between normal and cancer samples.

**Tool:** Metilene v0.2.6+galaxy0
**Input:** bedGraph files from all 6 samples
**Parameters:** q-value < 0.05

**Plot 1 — Mean Methylation Difference Distribution:**
![Metilene Plot 1](metilene_plot1_methylation_difference.png)

**Plot 2 — DMR Length Distribution (nucleotides):**
![Metilene Plot 2](metilene_plot2_dmr_length_nt.png)

**Plot 3 — DMR Length Distribution (CpG count):**
![Metilene Plot 3](metilene_plot3_dmr_length_cpg.png)

**Plot 4 — Q-value vs Mean Methylation Difference:**
![Metilene Plot 4](metilene_plot4_qvalue.png)

**Plot 5 — Group 1 vs Group 2 Mean Methylation:**
![Metilene Plot 5](metilene_plot5_group_comparison.png)

**Plot 6 — DMR Length (nt) vs DMR Length (CpGs):**
![Metilene Plot 6](metilene_plot6_dmr_scatter.png)

**Key findings:**
- Multiple DMRs identified between normal and cancer samples
- Both hypermethylated and hypomethylated regions found in cancer
- DMRs correspond to genes involved in cell identity and tumor suppression

---

##  Repository Structure
wgbs-breast-cancer-methylation/
│
├── Galaxy12-[MethylDackel on dataset 7_ fraction CpG].bedgraph

├── Galaxy13-[CpGIslands.bed].bed

├── Galaxy16-[heatmap2 on dataset 13].pdf

├── Galaxy50-[plotProfile on data 49_ Image].png
│
├── falco_dataset1_basic_stats.png

├──  falco_dataset1_per_base_quality.png

├──  falco_dataset1_per_seq_quality.png

├──  falco_dataset1_seq_content.png

├──  falco_dataset1_gc_content.png

├──  falco_dataset1_n_content.png

├──  falco_dataset1_seq_length.png

├── falco_dataset1_seq_duplication.png

├──  falco_dataset1_adapter_content.png

│
├── falco_dataset2_basic_stats.png
│
├──  igv_visualization.png

├──  methyldackel_cpg.png
│
├──  metilene_plot1_methylation_difference.png

├──  metilene_plot2_dmr_length_nt.png

├──  metilene_plot3_dmr_length_cpg.png

├──  metilene_plot4_qvalue.png

├──  metilene_plot5_group_comparison.png

├──  metilene_plot6_dmr_scatter.png
│
├──  plotprofile_single_sample.png

├──  plotprofile_all_samples.png
│
└──  README.md

---

##  How to Reproduce

1. Go to **usegalaxy.org** and create a free account
2. Create a new history named "WGBS"
3. Upload raw data from Zenodo:
https://zenodo.org/record/557099/files/subset_1.fastq
https://zenodo.org/record/557099/files/subset_2.fastq
4. Run **Falco** on both files
5. Run **bwameth** with hg38 reference (paired-end)
6. Run **MethylDackel** to extract CpG methylation fractions
7. Convert bedGraphs to bigWig
8. Run **computeMatrix** + **plotProfile** around CpG islands
9. Run **Metilene** to find DMRs
10. Full tutorial: https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html

---

## References

1. Lin et al. (2015) PLOS ONE. https://doi.org/10.1371/journal.pone.0118453
2. Galaxy WGBS Tutorial: https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html
3. Falco: https://falco.readthedocs.io
4. bwameth: https://github.com/brentp/bwa-meth
5. MethylDackel: https://github.com/dpryan79/MethylDackel
6. Metilene: https://www.bioinf.uni-leipzig.de/Software/metilene/
7. deepTools: https://deeptools.readthedocs.io
8. Galaxy: https://usegalaxy.org
