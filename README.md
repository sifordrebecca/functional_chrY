# Y Chromosome Callable Region and QV Analyses

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
