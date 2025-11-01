# SARS-CoV-2 Analysis Workflow - Complete Summary

## 📋 Quick Reference Card

### Project Overview
- **Project Name**: SARS-CoV-2 Genomic Surveillance Pipeline
- **Pipeline**: nf-core/viralrecon
- **Data Source**: CDC Benchmark Datasets
- **Samples**: 5 SARS-CoV-2 genomes
- **Platform**: Illumina paired-end sequencing
- **Protocol**: ARTIC amplicon

---

## 🔄 Complete Workflow (Step-by-Step)

### Phase 1: Setup (30 minutes)

```bash
# 1. Create project directory
mkdir sars-cov2-analysis
cd sars-cov2-analysis

# 2. Create conda environment
conda create --name nextflow nextflow
conda activate nextflow

# 3. Verify installation
nextflow help
```

### Phase 2: Data Acquisition (30-60 minutes)

```bash
# 4. Create sample list
cat > samples.txt << EOF
ERR5556343
SRR13500958
ERR5743893
ERR5181310
ERR5405022
EOF

# 5. Download data from SRA
for i in $(cat samples.txt); do 
    fastq-dump --split-files $i
done

# 6. Compress files
gzip *.fastq

# 7. Organize data
mkdir data
mv *.fastq.gz data/
```

### Phase 3: Samplesheet Preparation (5 minutes)

```bash
# 8. Download samplesheet script
wget -L https://raw.githubusercontent.com/nf-core/viralrecon/master/bin/fastq_dir_to_samplesheet.py

# 9. Generate samplesheet
python3 fastq_dir_to_samplesheet.py \
    data \
    samplesheet.csv \
    -r1 _1.fastq.gz \
    -r2 _2.fastq.gz

# 10. Verify samplesheet
cat samplesheet.csv
```

### Phase 4: Pipeline Execution (1-2 hours)

```bash
# 11. Ensure Docker is running
# Open Docker Desktop application

# 12. Run nf-core/viralrecon
nextflow run nf-core/viralrecon \
    -profile docker \
    --max_memory '12.GB' \
    --max_cpus 4 \
    --input samplesheet.csv \
    --outdir results/viralrecon \
    --protocol amplicon \
    --genome 'MN908947.3' \
    --primer_set artic \
    --primer_set_version 3 \
    --skip_kraken2 \
    --skip_assembly \
    --skip_pangolin \
    --skip_nextclade \
    --skip_asciigenome \
    --platform illumina \
    -resume
```

### Phase 5: Results Review (15 minutes)

```bash
# 13. Navigate to results
cd results/viralrecon

# 14. List outputs
ls

# 15. Open MultiQC report
open multiqc/multiqc_report.html

# 16. Check work directory size
cd ../..
du -sh work

# 17. Clean up (optional)
rm -rf work
```

### Phase 6: R Analysis (45 minutes)

```r
# 18. Open RStudio
# File → New Project → Existing Directory → sars-cov2-analysis

# 19. Install packages
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("RColorBrewer")

# 20. Load libraries
library(tidyverse)
library(ggplot2)

# 21. Load data
var <- read.csv("results/viralrecon/variants_long_table/combined_variants.csv")
var_tb <- as_tibble(var)

# 22. Explore data
var_tb %>% count(SAMPLE, sort = TRUE)
var_tb %>% group_by(SAMPLE) %>% summarize(mean(DP))

# 23. Create visualizations
ggplot(var_tb, aes(x = SAMPLE, y = DP, fill = SAMPLE)) + 
  geom_boxplot() + 
  ylim(0, 10000) +
  scale_fill_brewer(palette = "RdYlBu")

# 24. Save plots
ggsave("figures/depth_distribution.png", width = 10, height = 6, dpi = 300)

# Or run complete analysis script
source("variant_analysis.R")
```

---

## 📊 Key Results

### Quality Metrics Achieved

| Metric | Value | Status |
|--------|-------|--------|
| Samples Analyzed | 5 | ✅ Complete |
| Mean Coverage | 604-5,790× | ✅ Excellent |
| Genome Coverage (>10×) | 98-100% | ✅ Excellent |
| Mapping Rate | 94.59-100% | ✅ Excellent |
| Total SNPs | 20-32 per sample | ✅ Expected |
| Total INDELs | 1-6 per sample | ✅ Expected |
| Ambiguous Bases | <2.5% | ✅ High Quality |

### Output Files Generated

**Nextflow Pipeline:**
- ✅ 5 consensus genome sequences (FASTA)
- ✅ 5 variant call files (VCF/TSV)
- ✅ 5 alignment files (BAM)
- ✅ MultiQC quality report (HTML)
- ✅ Coverage statistics (TXT)
- ✅ Combined variants table (CSV)

**R Analysis:**
- ✅ 8 visualization plots (PNG)
- ✅ 4 summary tables (CSV)
- ✅ Depth statistics per sample
- ✅ Effect summary by sample
- ✅ Gene variant summary

---

## ⏱️ Time Estimates

| Phase | Task | Duration |
|-------|------|----------|
| 1 | Environment setup | 30 min |
| 2 | Data download | 30-60 min |
| 3 | Samplesheet creation | 5 min |
| 4 | Pipeline execution | 60-120 min |
| 5 | Results review | 15 min |
| 6 | R analysis | 45 min |
| **Total** | **Complete workflow** | **3-5 hours** |

*Note: Times vary based on internet speed and system specifications*

---

## 💾 Storage Requirements

| Component | Size | Notes |
|-----------|------|-------|
| Raw FASTQ files | ~5-10 GB | Compressed |
| Pipeline outputs | ~10-15 GB | Can be reduced |
| Work directory | ~4-5 GB | Can be deleted |
| R figures | ~20 MB | High resolution |
| Final project | ~15-25 GB | With work/ deleted |

---

## 🛠️ Software Requirements

### Core Tools
- **Nextflow**: ≥23.04.0
- **Java**: ≥11
- **Docker**: Latest stable
- **Conda**: ≥4.12.0
- **SRA Toolkit**: ≥3.0.0

### R Packages
- tidyverse (includes ggplot2, dplyr)
- ggpubr
- RColorBrewer

### System Specs
- **RAM**: 8GB minimum (12GB recommended)
- **CPU**: 4 cores minimum
- **Storage**: 50GB free space
- **OS**: Linux, macOS, or Windows (with WSL)

---

## 📈 Skills Demonstrated

### Bioinformatics
- ✅ Viral genome assembly
- ✅ Quality control assessment
- ✅ Variant calling and annotation
- ✅ Coverage analysis
- ✅ Reference-based alignment

### Computational
- ✅ Workflow management (Nextflow)
- ✅ Containerization (Docker)
- ✅ Shell scripting (Bash)
- ✅ Statistical analysis (R)
- ✅ Data visualization (ggplot2)

### Research
- ✅ Public data retrieval
- ✅ Reproducible research
- ✅ Documentation
- ✅ Version control (Git)
- ✅ Scientific communication

---

## 🎯 Learning Outcomes

By completing this project, you will:

1. **Master nf-core pipelines** - Industry-standard workflows
2. **Understand viral genomics** - SARS-CoV-2 analysis
3. **Work with public data** - SRA/ENA databases
4. **Develop bioinformatics skills** - QC, alignment, variant calling
5. **Create reproducible research** - Documentation and version control
6. **Visualize genomic data** - R and ggplot2
7. **Apply to public health** - Surveillance applications

---

## 📚 Key Commands Cheatsheet

```bash
# Environment
conda activate nextflow
conda deactivate

# Pipeline
nextflow run nf-core/viralrecon -profile docker [params]
nextflow log
nextflow clean -f

# Data
fastq-dump --split-files SRR13500958
gzip *.fastq
du -sh directory/

# Results
ls -lh results/viralrecon/
cat samplesheet.csv
open multiqc/multiqc_report.html

# Git
git add .
git commit -m "message"
git push origin main
```

```r
# R Basics
var <- read.csv("file.csv")
var_tb <- as_tibble(var)
head(var_tb)
summary(var_tb)

# Data Manipulation
var_tb %>% count(SAMPLE)
var_tb %>% filter(DP >= 500)
var_tb %>% group_by(SAMPLE) %>% summarize(mean(DP))

# Visualization
ggplot(var_tb, aes(x=SAMPLE, y=DP)) + geom_boxplot()
ggsave("plot.png", width=10, height=6, dpi=300)
```

---

## 🚨 Common Issues & Quick Fixes

| Problem | Solution |
|---------|----------|
| `nextflow: command not found` | `conda activate nextflow` |
| `No space left on device` | `rm -rf work/` |
| `Java heap space error` | `export NXF_OPTS="-Xmx4G"` |
| `Docker not running` | Open Docker Desktop |
| `Cannot find samplesheet` | Check file path and name |
| `R package not found` | `install.packages("package")` |

---

## ✅ Success Checklist

**Pipeline Completion:**
- [ ] All 5 samples processed
- [ ] MultiQC report generated
- [ ] Consensus sequences created
- [ ] No error messages in log
- [ ] Coverage >98% for all samples

**R Analysis Completion:**
- [ ] Variant data loaded successfully
- [ ] 8 plots generated in figures/
- [ ] 4 tables exported to tables/
- [ ] No errors in script execution
- [ ] Plots are publication-quality

**Documentation:**
- [ ] README completed with your info
- [ ] Images added to repository
- [ ] All scripts are executable
- [ ] Repository is public on GitHub
- [ ] .gitignore is working properly

---

## 📝 Citation

**This Project:**
```
[Your Name]. (2024). SARS-CoV-2 Genomic Surveillance Pipeline. 
GitHub repository: https://github.com/yourusername/sars-cov2-analysis
```

**CDC Datasets:**
```
Li X, Hagey JV, Moran EJ, et al. (2022). Benchmark datasets for 
SARS-CoV-2 surveillance bioinformatics. PeerJ 10:e13821.
https://doi.org/10.7717/peerj.13821
```

**nf-core/viralrecon:**
```
Ewels PA, et al. (2020). The nf-core framework for community-curated 
bioinformatics pipelines. Nature Biotechnology 38:276-278.
https://doi.org/10.1038/s41587-020-0439-x
```

