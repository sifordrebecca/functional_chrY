# ============================================================
# EULER / VENN DIAGRAMS FOR CALLABLE-MASK OVERLAP
# ============================================================
#Written by Siford / Finalized May 2026 / sifordrebecca@gmail.com

#Note the Euler / venn diagrams are very finicky - 
#check PDF permissions if it fails more than once if using a Mac 

#Each Euler diagram shows, within one CHM13 Y sequence class, how many callable bases are:
#unique to each mask,shared by two masks, or shared by all three - also output to a chart for easier viewing 

#Blue = Pangenome, Orange = Short-read T2T, Purple = Short-read GRCh37

#note that the overlap between the 3 masks was accomplished with bedtools multiinter



library(eulerr)
library(tidyverse)
# ------------------------------------------------------------
# Read bedtools multiinter output
# ------------------------------------------------------------

multiinter_dir <- "~/Downloads/HENG_HS1_POZ_overlap"

read_multiinter <- function(file) {
  
  x <- read.table(
    file.path(multiinter_dir, file),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(x) <- c(
    "chr",
    "start",
    "end",
    "count",
    "members",
    "pm151b",
    "hs1",
    "poznik"
  )
  
  x %>%
    mutate(
      bp = end - start,
      members = as.character(members)
    )
}

# ------------------------------------------------------------
# Read individual sequence classes
# ------------------------------------------------------------

m_xdr    <- read_multiinter("X-DEG.multiinter.tsv")
m_other  <- read_multiinter("OTHER.multiinter.tsv")
m_ampl   <- read_multiinter("AMPL.multiinter.tsv")
m_xtr1   <- read_multiinter("XTR1.multiinter.tsv")
m_xtr2   <- read_multiinter("XTR2.multiinter.tsv")
m_sat    <- read_multiinter("SAT.multiinter.tsv")
m_dyz17  <- read_multiinter("DYZ17.multiinter.tsv")
m_dyz19  <- read_multiinter("DYZ19.multiinter.tsv")
m_cen    <- read_multiinter("CEN.multiinter.tsv")
m_het    <- read_multiinter("HET.multiinter.tsv")

# Combine classes to match callable heatmap (Figure 5h)
m_xtr <- bind_rows(
  m_xtr1,
  m_xtr2
)

m_sat_combined <- bind_rows(
  m_sat,
  m_dyz17
)

# ------------------------------------------------------------
# Convert multiinter output into Euler counts
# ------------------------------------------------------------

get_euler_counts <- function(dat) {
  
  totals <- dat %>%
    group_by(members) %>%
    summarise(
      bp = sum(bp),
      .groups = "drop"
    )
  
  get_bp <- function(code) {
    
    value <- totals$bp[
      totals$members == code
    ]
    
    if (length(value) == 0)
      return(0)
    
    sum(value)
  }
  
  c(
    "Pangenome" =
      get_bp("1"),
    
    "Short-read_T2T" =
      get_bp("2"),
    
    "Short-read_GRCh37" =
      get_bp("3"),
    
    "Pangenome&Short-read_T2T" =
      get_bp("1,2"),
    
    "Pangenome&Short-read_GRCh37" =
      get_bp("1,3"),
    
    "Short-read_T2T&Short-read_GRCh37" =
      get_bp("2,3"),
    
    "Pangenome&Short-read_T2T&Short-read_GRCh37" =
      get_bp("1,2,3")
  )
}

# ------------------------------------------------------------
# Counts for each final sequence class
# ------------------------------------------------------------

counts_xdr <- get_euler_counts(m_xdr)

counts_other <- get_euler_counts(m_other)

counts_ampl <- get_euler_counts(m_ampl)

counts_xtr <- get_euler_counts(m_xtr)

counts_sat <- get_euler_counts(m_sat_combined)

counts_dyz19 <- get_euler_counts(m_dyz19)

counts_cen <- get_euler_counts(m_cen)

counts_het <- get_euler_counts(m_het)

# ------------------------------------------------------------
# Put all counts in one list
# ------------------------------------------------------------

euler_counts <- list(
  XDR = counts_xdr,
  OTHER = counts_other,
  AMPL = counts_ampl,
  XTR = counts_xtr,
  SAT = counts_sat,
  DYZ19 = counts_dyz19,
  `CEN-DYZ3` = counts_cen,
  HET = counts_het
)

# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------

set_colors <- c(
  "Pangenome" = "#5DA5DA",
  "Short-read_T2T" = "#E69F00",
  "Short-read_GRCh37" = "#C060C0"
)

# ------------------------------------------------------------
# Plotting function
# ------------------------------------------------------------

plot_euler_class <- function(counts, title) {
  
  # Remove overlap categories that are exactly zero
  counts_use <- counts[counts > 0]
  
  fit <- euler(counts_use)
  
  # Determine which individual sets actually occur
  present_sets <- unique(
    unlist(
      strsplit(
        names(counts_use),
        "&",
        fixed = TRUE
      )
    )
  )
  
  present_colors <- set_colors[present_sets]
  
  plot(
    fit,
    
    fills = list(
      fill = present_colors,
      alpha = 0.65
    ),
    
    edges = list(
      col = "black",
      lwd = 2.2
    ),
    
    labels = list(
      fontsize = 11
    ),
    
    quantities = list(
      type = "counts",
      fontsize = 10
    ),
    
    main = title
  )
}

# ============================================================
# Plot each sequence class - edit in Illustrator 
# ============================================================

plot_euler_class(
  counts_xdr,
  "XDR"
)

plot_euler_class(
  counts_other,
  "OTHER"
)

plot_euler_class(
  counts_ampl,
  "AMPL"
)

plot_euler_class(
  counts_xtr,
  "XTR"
)

plot_euler_class(
  counts_sat,
  "SAT"
)

plot_euler_class(
  counts_dyz19,
  "DYZ19"
)

plot_euler_class(
  counts_cen,
  "CEN-DYZ3"
)

plot_euler_class(
  counts_het,
  "HET"
)

# ============================================================
# Save overlap counts as a table
# ============================================================


euler_summary <- bind_rows(
  lapply(
    names(euler_counts),
    function(cat) {
      
      tibble(
        CATEGORY = cat,
        overlap = names(euler_counts[[cat]]),
        bp = as.numeric(euler_counts[[cat]])
      )
    }
  )
)

write.table(
  euler_summary,
  "callable_mask_Euler_overlap_counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Optional Mb version
euler_summary_mb <- euler_summary %>%
  mutate(
    Mb = bp / 1e6
  )

write.table(
  euler_summary_mb,
  "callable_mask_Euler_overlap_counts_Mb.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
