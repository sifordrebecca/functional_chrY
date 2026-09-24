#TURTH ASSEMBLY CN VS NORMALIZED TSPY1 READ DEPTH - final + linear model for prediction 

# Files
master_file <- "TSPY1_depth_normalized_master.withGenome.tsv"
truth_file  <- "tspy_truthCN.csv"

# Read
master <- read.table(master_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
truth  <- read.csv(truth_file, stringsAsFactors=FALSE)

# Prep truth
truth2 <- truth[, c("sample", "TSPY")]
names(truth2) <- c("sample", "truthCN")

# Merge
df <- merge(master, truth2, by="sample", all=FALSE)

# Numeric cleanup
cols_to_num <- c("truthCN","tspy1_norm_y","tspy1_per_Mreads","tspy1_norm_genome")
for (cc in cols_to_num) df[[cc]] <- as.numeric(df[[cc]])
df <- df[is.finite(df$truthCN), ]

# Fit models (CN on Y)
m_y     <- lm(truthCN ~ tspy1_norm_y, data=df)
m_mread <- lm(truthCN ~ tspy1_per_Mreads, data=df)
m_gen   <- lm(truthCN ~ tspy1_norm_genome, data=df)

r2_y     <- summary(m_y)$r.squared
r2_mread <- summary(m_mread)$r.squared
r2_gen   <- summary(m_gen)$r.squared

cat("R^2 (DDX3Y ratio):", round(r2_y, 4), "\n")
cat("R^2 (per M reads):", round(r2_mread, 4), "\n")
cat("R^2 (genome depth):", round(r2_gen, 4), "\n")

# Plots (3 separate)
plot(df$tspy1_norm_y, df$truthCN,
     xlab="tspy1_norm_y (TSPY1/DDX3Y)",
     ylab="Truth TSPY copy number",
     main=paste0("CN vs tspy1_norm_y (R²=", round(r2_y, 3), ")"),
     pch=16)
abline(m_y, lwd=2)

plot(df$tspy1_per_Mreads, df$truthCN,
     xlab="tspy1_per_Mreads (depth per million mapped reads)",
     ylab="Truth TSPY copy number",
     main=paste0("CN vs tspy1_per_Mreads (R²=", round(r2_mread, 3), ")"),
     pch=16)
abline(m_mread, lwd=2)

plot(df$tspy1_norm_genome, df$truthCN,
     xlab="tspy1_norm_genome (TSPY1 / autosome mean depth)",
     ylab="Truth TSPY copy number",
     main=paste0("CN vs tspy1_norm_genome (R²=", round(r2_gen, 3), ")"),
     pch=16)
abline(m_gen, lwd=2)


#make chart 
# Summary table: normalization method + R² + slope + intercept + slope p-value

master_file <- "TSPY1_depth_normalized_master.withGenome.tsv"
truth_file  <- "tspy_truthCN.csv"

master <- read.table(master_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
truth  <- read.csv(truth_file, stringsAsFactors=FALSE)

truth2 <- truth[, c("sample", "TSPY")]
names(truth2) <- c("sample", "truthCN")

df <- merge(master, truth2, by="sample", all=FALSE)

# numeric cleanup
num_cols <- c("truthCN","tspy1_norm_y","tspy1_per_Mreads","tspy1_norm_genome")
for (cc in num_cols) df[[cc]] <- as.numeric(df[[cc]])
df <- df[is.finite(df$truthCN), ]

make_row <- function(method_name, predictor_col) {
  d <- df[is.finite(df[[predictor_col]]), ]
  fit <- lm(truthCN ~ d[[predictor_col]], data=d)
  s <- summary(fit)
  
  data.frame(
    method = method_name,
    n = nrow(d),
    intercept = unname(coef(fit)[1]),
    slope = unname(coef(fit)[2]),
    r_squared = s$r.squared,
    slope_p_value = s$coefficients[2, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}

res <- rbind(
  make_row("Y baseline (TSPY1/DDX3Y)", "tspy1_norm_y"),
  make_row("Library size (per million reads)", "tspy1_per_Mreads"),
  make_row("Genome depth (TSPY1/autosome mean)", "tspy1_norm_genome")
)

print(res)
# Optional: nicer formatting
res$r_squared <- round(res$r_squared, 4)
res$intercept <- round(res$intercept, 4)
res$slope <- round(res$slope, 4)
res$slope_p_value <- signif(res$slope_p_value, 3)
print(res)

#----Linear Model write up below 
# using autosomal depth as choice for now to convert depth to CN for gtex samples 
#fitting the model with CN as the response and depth metric as the predictor // extract the slope (m) and intercept (b) and print them as a formula

# Fit calibration: CN = m * depth_metric + b
fit <- lm(truthCN ~ tspy1_norm_genome, data = df)   # replace depth_metric with column name

# Extract coefficients
b <- coef(fit)[1]   # intercept
m <- coef(fit)[2]   # slope

# Print the formula - for autosomal norm depth for now 
cat(sprintf("CN_hat = %.4f * normalized depth + %.4f\n", m, b))

# Define a reusable function (so you can apply it to GTEx later)
predict_CN <- function(depth_value) {
  m * depth_value + b
}

#will use normalized depth from autosomes over tspy1 from each gtex sample to get predicted CN using formula 




