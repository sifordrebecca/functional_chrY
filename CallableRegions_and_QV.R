#####    Final T2T Y CALLABLE REGIONS:     ####
#####   1. Callable fractions (Figure 5h)  ####
#####   2. Mean QV (Figure 5i)             ####

###############################################
###############################################

#Written by Siford / Finalized Sep 24.2026 / sifordrebecca@gmail.com

###############################################
###############################################
# 1. Callable Fractions across Y sequence classes and Y masks # (+Venn Diagrams)

#this script section:  
#uses the CHM13 sequence-class annotation: https://github.com/marbl/CHM13
#intersects it with the three global masks: https://zenodo.org/records/14911560, https://pmc.ncbi.nlm.nih.gov/articles/PMC4032117/#SM, https://www.illumina.com/science/genomics-research/articles/identifying-genomic-regions-with-high-quality-single-nucleotide-.html
#combines DYZ17 + SAT → SAT, combines XTR1 + XTR2 → XTR, renames X-DEG → XDR, renames CEN → CEN-DYZ3
#shows callable fraction by color (blue scale)
#shows retained Mb inside each tile
#makes true 0 / missing combinations grey
#puts total mask bp available (Mb) in a separate row at the bottom 

library(tidyverse)
library(GenomicRanges)
library(IRanges)
library(scales)
library(patchwork)

# ============================================================
# Files: global mask bed files + chm13 Y 
# ============================================================

annotation_file <- "chm13v2.0_chrY_sequence_class_v1.bed"
#note that the Y has been excised out of the XY major bed file before loading into R 

#excluded masks produced using hg38, pm151a, see methods 
mask_files <- c(
  "Short read - GRCh37" = "Poznik_2013_S1b.chrY.callable.chm13.merged.bed",
  "Short read - T2T"    = "hs1.Y_mask.bed",
  "Pangenome"           = "pm151b_v3.easy.chrY.bed" #more stringent than a 
)

# ============================================================
# Functions
# ============================================================

read_bed_gr <- function(file) {
  
  x <- read.table(
    file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  x <- x[, 1:3]
  colnames(x) <- c("chr", "start", "end")
  
  GRanges(
    seqnames = x$chr,
    ranges = IRanges(
      start = x$start + 1,
      end = x$end
    )
  ) %>%
    reduce()
}


format_mb_label <- function(x_bp) {
  
  x_mb <- x_bp / 1e6
  
  sapply(x_mb, function(x) {
    
    if (is.na(x)) return("")
    
    if (x == 0) return("0")
    
    if (x < 0.001)
      return(formatC(
        x,
        format = "e",
        digits = 1
      ))
    
    sprintf("%.3f", x)
  })
}

# ============================================================
# Read CHM13 sequence-class annotation
# ============================================================

annot <- read.table(
  annotation_file,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

annot <- annot[, 1:4]

colnames(annot) <- c(
  "chr",
  "start",
  "end",
  "CATEGORY"
)

# Harmonize sequence classes / based on current representation of Y chr 

annot <- annot %>%
  mutate(
    CATEGORY = case_when(
      CATEGORY == "X-DEG" ~ "XDR",
      CATEGORY == "CEN" ~ "CEN-DYZ3",
      CATEGORY %in% c("DYZ17", "SAT") ~ "SAT",
      CATEGORY %in% c("XTR1", "XTR2") ~ "XTR",
      TRUE ~ CATEGORY
    )
  )

# Final classes shown in plot

keep_cats <- c(
  "XDR",
  "OTHER",
  "AMPL",
  "XTR",
  "SAT",
  "DYZ19",
  "CEN-DYZ3",
  "HET"
)

annot <- annot %>%
  filter(CATEGORY %in% keep_cats)

# Convert annotation to GRanges

annot_gr <- GRanges(
  seqnames = annot$chr,
  ranges = IRanges(
    start = annot$start + 1,
    end = annot$end
  ),
  CATEGORY = annot$CATEGORY
)

# Reduce separately within each sequence class

class_gr <- split(
  annot_gr,
  annot_gr$CATEGORY
)

class_gr <- lapply(
  class_gr,
  reduce
)

# ============================================================
# Read global masks
# ============================================================

mask_gr <- lapply(
  mask_files,
  read_bed_gr
)

# ============================================================
# Calculate callable fraction by sequence class
# ============================================================

callable <- expand_grid(
  CATEGORY = keep_cats,
  MASK = names(mask_gr)
) %>%
  
  rowwise() %>%
  
  mutate(
    
    bp_total = {
      gr <- class_gr[[CATEGORY]]
      sum(width(gr))
    },
    
    bp_in_easy = {
      gr_class <- class_gr[[CATEGORY]]
      gr_mask  <- mask_gr[[MASK]]
      
      ov <- intersect(
        gr_class,
        gr_mask
      )
      
      sum(width(ov))
    },
    
    frac_in_easy = bp_in_easy / bp_total
    
  ) %>%
  
  ungroup()

# ============================================================
# Prep to plot  
# ============================================================

callable <- callable %>%
  mutate(
    
    CATEGORY = factor(
      CATEGORY,
      levels = rev(keep_cats)
    ),
    
    MASK = factor(
      MASK,
      levels = c(
        "Short read - GRCh37",
        "Short read - T2T",
        "Pangenome"
      )
    ),
    
    # Keep 0 values grey instead of putting them on the blue gradient
    fraction_plot = case_when(
      is.na(frac_in_easy) ~ NA_real_,
      frac_in_easy == 0 ~ NA_real_,
      TRUE ~ frac_in_easy
    ),
    
    mb_label = format_mb_label(
      bp_in_easy
    ),
    
    text_col = case_when(
      !is.na(frac_in_easy) &
        frac_in_easy > 0.60 ~ "white",
      TRUE ~ "black"
    )
  )

# ============================================================
# Main callable-fraction heatmap
# ============================================================

p_main <- ggplot(
  callable,
  aes(
    x = MASK,
    y = CATEGORY
  )
) +
  
  geom_tile(
    aes(fill = fraction_plot),
    width = 1,
    height = 1,
    color = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      label = mb_label,
      color = text_col
    ),
    size = 4
  ) +
  
  scale_color_identity() +
  
  scale_fill_gradientn(
    colors = c(
      "#f7fbff",
      "#6baed6",
      "#2171b5",
      "#08306b"
    ),
    limits = c(0, 1),
    oob = scales::squish,
    na.value = "grey80",
    name = "Callable\nfraction"
  ) +
  
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Callable Fraction of Sequence by Mask and\nY Chromosome Sequence Class",
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    
    legend.position = "right"
  )

# ============================================================
# Total mask bp available
# ============================================================

mask_totals <- tibble(
  MASK = names(mask_gr),
  total_bp = sapply(
    mask_gr,
    function(x) sum(width(x))
  )
) %>%
  
  mutate(
    
    MASK = factor(
      MASK,
      levels = c(
        "Short read - GRCh37",
        "Short read - T2T",
        "Pangenome"
      )
    ),
    
    label = sprintf(
      "%.3f",
      total_bp / 1e6
    )
  )

# ============================================================
# Bottom total-mask row
# ============================================================

p_bottom <- ggplot(
  mask_totals,
  aes(
    x = MASK,
    y = 1
  )
) +
  
  geom_tile(
    fill = "grey90",
    width = 1,
    height = 1,
    color = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(label = label),
    size = 4
  ) +
  
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = 1,
    labels = "Total mask bp\navailable (Mb)",
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# ============================================================
# Combine
# ============================================================

callable_figure <- p_main / p_bottom +
  plot_layout(
    heights = c(8, 1.25)
  )

callable_figure

# ============================================================
# Save
# ============================================================

ggsave(
  "callable_fraction_by_mask_sequence_class.pdf",
  callable_figure,
  width = 6,
  height = 6.5
)

ggsave(
  "callable_fraction_by_mask_sequence_class.png",
  callable_figure,
  width = 6,
  height = 6.5,
  dpi = 300
)

###############################################
#Note the Euler / venn diagrams are very finicky - 
#check PDF permissions if it fails more than once if using a Mac 

#see different script for this withing same repository 


###############################################
###############################################
# 2. Mean QV across Y sequence classes and Y masks, and mean QV for all regions by mask

# This script relies on a ERRBASE masked summary file, code to configure this file can be found here https://github.com/nhansen/chrY_callable_analysis
#Remember to EXCLUDE HG002 from your data before making ERRBASE file

library(tidyverse)
library(scales)
library(patchwork)

# ============================================================
# 1. Read ERRBASE-masked summary file
# ============================================================

qv <- read.table(
  "summary_category_stats_minalign1000_refannot_errorbase_masking_no_HG002.txt",
  header = FALSE,
  stringsAsFactors = FALSE,
  na.strings = "NA"
)

colnames(qv) <- c(
  "CATEGORY",
  "MASK",
  "TOTAL_BASES",
  "ALIGNED_BASES",
  "NUM_DISCREPANCIES",
  "ERROR_RATE",
  "PCT_ALIGNED",
  "QV"
)

# ============================================================
# 2. Sequence classes and masks
# ============================================================

keep_cats <- c(
  "XDR",
  "OTHER",
  "AMPL",
  "XTR",
  "SAT",
  "DYZ19",
  "CEN-DYZ3",
  "HET"
)

mask_names <- c(
  "poznik2013" = "Short read - GRCh37",
  "hs1"        = "Short read - T2T",
  "pm151bv3"   = "Pangenome"
)

# ============================================================
# 3. Prep QV heatmap data
# ============================================================

qv_plot <- qv %>%
  filter(
    CATEGORY %in% keep_cats,
    MASK %in% names(mask_names)
  ) %>%
  mutate(
    MASK = recode(MASK, !!!mask_names),
    
    MASK = factor(
      MASK,
      levels = c(
        "Short read - GRCh37",
        "Short read - T2T",
        "Pangenome"
      )
    ),
    
    CATEGORY = factor(
      CATEGORY,
      levels = rev(keep_cats)
    ),
    
    # Inf = no observed error.
    # Display at 60 for the color scale, but label as 60+.
    # Finite values >60 remain unchanged.
    QV_fill = case_when(
      is.infinite(QV) ~ 60,
      TRUE ~ QV
    ),
    
    QV_label = case_when(
      is.na(QV) ~ "",
      is.infinite(QV) ~ "60+",
      TRUE ~ sprintf("%.1f", QV)
    ),
    
    text_col = ifelse(
      !is.na(QV_fill) & QV_fill > 53,
      "white",
      "black"
    )
  )

# ============================================================
# 4. White boxes
# ============================================================

#White boxes = Boxes are changed to white when not enough sequence data is availble to compute, the mask retained some sequence, but less than 1,000-bp minimum, so QV was not calculated.

# T2T x CEN-DYZ3, Pangenome x DYZ19 hardcoded 
#
# Other missing QVs remain grey.

white_boxes <- qv_plot %>%
  filter(
    (MASK == "Short read - T2T" & CATEGORY == "CEN-DYZ3") |
      (MASK == "Pangenome" & CATEGORY == "DYZ19")
  )

qv_main <- qv_plot %>%
  mutate(
    special_white =
      (MASK == "Short read - T2T" & CATEGORY == "CEN-DYZ3") |
      (MASK == "Pangenome" & CATEGORY == "DYZ19"),
    
    QV_fill = ifelse(
      special_white,
      NA_real_,
      QV_fill
    ),
    
    QV_label = ifelse(
      special_white,
      "",
      QV_label
    )
  )

# ============================================================
# 5. Main QV heatmap
# ============================================================
#Grey boxes = the mask retained no callable sequence for that sequence class, so there was nothing to evaluate.

p_main <- ggplot(
  qv_main,
  aes(x = MASK, y = CATEGORY)
) +
  
  # Main tiles
  # Missing values are grey
  geom_tile(
    aes(fill = QV_fill),
    width = 1,
    height = 1,
    color = "white",
    linewidth = 0.5
  ) +
  
  # Overlay the two intentionally white boxes
  geom_tile(
    data = white_boxes,
    fill = "white",
    width = 1,
    height = 1,
    color = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      label = QV_label,
      color = text_col
    ),
    size = 4
  ) +
  
  scale_color_identity() +
  
  scale_fill_gradientn(
    colors = c(
      "#f7fcf5",
      "#c7e9c0",
      "#74c476",
      "#238b45",
      "#00441b"
    ),
    limits = c(20, 60),
    oob = scales::squish,
    na.value = "grey80",
    name = "Mean QV"
  ) +
  
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Mean QV",
    subtitle =
      "QV = -10log10(summed discrepancies /\nsummed aligned bases); HG002 excluded",
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Remove little dashes after class names to save space 
    axis.ticks.y = element_blank(),
    
    axis.line = element_blank(),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11,
      lineheight = 1.1
    ),
    
    legend.position = "right"
  )

# ============================================================
# 6. CALCULATE overall QV for all regions together by mask 
# ============================================================
#
# QV_avg =
# -10 * log10(
#   total discrepancies / total aligned bases
# )
#
# This is NOT the arithmetic mean of the QV values, cant take the mean of logged values
# It pools discrepancies and aligned bases 
# ============================================================

overall <- qv %>%
  filter(
    MASK %in% c(
      "poznik2013",
      "hs1",
      "pm151bv3"
    )
  ) %>%
  
  group_by(MASK) %>%
  
  summarise(
    total_discrepancies =
      sum(NUM_DISCREPANCIES, na.rm = TRUE),
    
    total_aligned =
      sum(ALIGNED_BASES, na.rm = TRUE),
    
    QV = case_when(
      total_aligned == 0 ~ NA_real_,
      total_discrepancies == 0 ~ Inf,
      TRUE ~
        -10 * log10(
          total_discrepancies /
            total_aligned
        )
    ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    MASK = recode(
      MASK,
      "poznik2013" = "Short read - GRCh37",
      "hs1"        = "Short read - T2T",
      "pm151bv3"   = "Pangenome"
    ),
    
    MASK = factor(
      MASK,
      levels = c(
        "Short read - GRCh37",
        "Short read - T2T",
        "Pangenome"
      )
    ),
    
    label = case_when(
      is.infinite(QV) ~ "60+",
      TRUE ~ sprintf("%.2f", QV)
    )
  )

# ============================================================
# 7. Check calculated overall values
# ============================================================

print(overall)

# Expected approximately:
#
# Short read - GRCh37   49.30
# Short read - T2T      49.33
# Pangenome             64.91

# ============================================================
# 8. Grey "Mean QV all regions" bottom row
# ============================================================

p_bottom <- ggplot(
  overall,
  aes(x = MASK, y = 1)
) +
  
  geom_tile(
    fill = "grey90",
    width = 1,
    height = 1,
    color = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(label = label),
    color = "black",
    size = 4
  ) +
  
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = 1,
    labels = "Mean QV\nall regions",
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# ============================================================
# 9. Combine
# ============================================================

qv_figure <- p_main / p_bottom +
  plot_layout(
    heights = c(8, 1.25)
  )

qv_figure

# ============================================================
# 10. Save
# ============================================================

ggsave(
  "mean_QV_by_mask_sequence_class_ERRBASE_masked_nohg002.pdf",
  qv_figure,
  width = 5,
  height = 7
)

ggsave(
  "mean_QV_by_mask_sequence_class_ERRBASE_masked_nohg002.png",
  qv_figure,
  width = 5,
  height = 7,
  dpi = 300
)
