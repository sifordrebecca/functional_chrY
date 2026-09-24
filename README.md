# Y Chromosome Callable Region and QV Analyses and TSPY Analyses

This repository contains scripts used to compare three global Y chromosome callable-region masks and evaluate short-read variant accuracy across T2T-CHM13v2Y sequence classes.
The three masks compared are:
- **Short read - GRCh37**: Poznik et al. callable Y mask
- **Short read - T2T**: T2T-CHM13v2Y short-read-based mask
- **Pangenome**: pm151b pangenome-based mask from Heng Li 

Sequence classes are based on the T2T-CHM13v2Y annotation and are summarized as:
`XDR`, `OTHER`, `AMPL`, `XTR`, `SAT`, `DYZ19`, `CEN-DYZ3`, and `HET`.

---
## Callable Regions
Compare how much Y chromosome sequence is retained by each global callable-region mask
See 'CallableRegions_and_QV.R' and 'callable_mask_overlap_euler.R'

### Analysis
Each mask was intersected with the T2T-CHM13v2Y sequence-class annotation. For each sequence class, the analysis calculates:
- total annotated sequence
- callable bases retained by each mask
- fraction of the sequence class retained

For plotting, related annotation classes were combined:
- `X-DEG` → `XDR`
- `XTR1 + XTR2` → `XTR`
- `SAT + DYZ17` → `SAT`
- `CEN` → `CEN-DYZ3`

The heatmap displays callable fraction by color and retained sequence in Mb within each tile. The bottom row reports the total sequence available in each global mask.

---
## QV
Measure short-read variant accuracy within the regions retained by each mask.
See 'CallableRegions_and_QV.R'

### Analysis

Short-read variant calls were used to construct sample-specific Y chromosome pseudogenomes. These were compared with the corresponding Verkko2 Y chromosome assemblies.
Bases annotated as `ERRBASE` in the sample-specific Verkko2 assemblies were masked before comparison so that low-confidence assembly positions were not included in the QV calculation. HG002 was excluded.

QV is calculated from discrepancy rate:
`QV = -10 × log10(discrepancies / aligned bases)`

For mean QV, discrepancies and aligned bases are pooled before calculating QV rather than averaging individual QV values:
`QV_mean = -10 × log10(sum(discrepancies) / sum(aligned bases))`

Additional plotting:
- categories with insufficient aligned sequence are not assigned a QV
- `Inf` QV values indicate no observed discrepancies and are displayed as `60+`
- finite QV values greater than 60 are retained
- grey and white boxes distinguish unavailable and insufficiently supported comparisons


---
## Venn / Euler Diagrams
Show how callable bases are shared or unique among the three global masks within each Y chromosome sequence class.

### Analysis
Mask overlaps were generated with `bedtools multiinter`. Each interval is assigned to one of seven possible overlap categories:
- Pangenome only
- Short read - T2T only
- Short read - GRCh37 only
- Pangenome + T2T
- Pangenome + GRCh37
- T2T + GRCh37
- all three masks

Base-pair lengths are summed within each overlap category and visualized with Euler diagrams.
As in the callable-region analysis:

- `XTR1 + XTR2` are combined as `XTR`
- `SAT + DYZ17` are combined as `SAT`

These diagrams summarize the genomic overlap among the global masks and are independent of the sample-specific QV analysis.


--
# TSPY1 Short-Read Depth and Copy Number Calibration 

This repository contains two connected workflows for estimating TSPY copy number from short-read whole-genome sequencing data. The first workflow calculates TSPY1 read depth and multiple normalization metrics from BAM/CRAM files. The second workflow compares those normalized depth values with assembly-derived TSPY copy number and fits a linear calibration model that can be used to estimate copy number in new samples, in this case from GTEX.

The workflow was developed for GRCh38-aligned whole-genome sequencing data.

1. Calculate short-read depth across the **TSPY1** locus and several normalization regions from BAM/CRAM files.
2. Compare normalized TSPY1 read depth with assembly-derived TSPY copy number.
3. Fit a linear calibration model that can be applied to additional samples, including GTEx.

---

## Workflow 1: Calculate TSPY1 Read Depth and Normalization Metrics

The first workflow is a SLURM array job that processes one BAM/CRAM file per task.
See `gtex_generate_depthmetrics.sh`.

For each sample, mean read depth is calculated across the GRCh38 TSPY1 locus using `samtools depth`.
The workflow also calculates several normalization baselines:
- DDX3Y mean depth
- autosomal mean depth across chromosomes 1–22
- total mapped reads (excluded)
- mean depth across X-degenerate Y chromosome regions

The main normalization metrics are:
- `tspy1_norm_y` = TSPY1 depth / DDX3Y depth
- `tspy1_norm_genome` = TSPY1 depth / autosomal mean depth
- `tspy1_norm_xdr` = TSPY1 depth / XDR mean depth

### Required inputs

The script requires:
- a list of BAM/CRAM paths, one file per line
- GRCh38 reference FASTA; the reference must match the reference used to generate the alignments
- TSPY1 BED file
- DDX3Y BED file
- autosomal BED file
- XDR BED file

---

## Workflow 2: Calibrate Normalized TSPY1 Depth to TSPY Copy Number

The second workflow is an R analysis that compares normalized TSPY1 depth with assembly-derived TSPY copy number.
See `CopyNumber_vs_Depth.R`.

The goal is to determine how well short-read depth predicts assembly-based copy number and to generate a linear equation that can be applied to GTEx samples.

The analysis combines:
- assembly-derived TSPY copy-number truth
- normalized TSPY1 short-read depth

Several normalization approaches are compared. TSPY1 depth normalized by autosomal mean depth was selected as the primary model: `TSPY copy number ~ TSPY1 / autosomal mean depth`

The script generates regression statistics and scatterplots and fits the calibration model used to estimate TSPY copy number in additional short-read datasets.

---

## Overall Workflow

Workflow 1 generates normalized TSPY1 depth measurements.
Workflow 2 uses those measurements to build and evaluate the depth-to-copy-number calibration model.

```text
BAM/CRAM
   ↓
TSPY1 read depth
   ↓
Normalized TSPY1 depth
   ↓
Comparison with assembly-derived copy number
   ↓
Linear calibration model
   ↓
Estimated TSPY copy number in new samples
