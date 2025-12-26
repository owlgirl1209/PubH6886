# Simulate per-sample cell-level CSVs where a few random (marker, cell-type)
# pairs have marker intensities that are positively correlated with subject survival.
#
# Output:
#  - {SampleID}_cells.csv files (one per sample)
#  - survival_associations.csv listing which (marker, cell-type) pairs are correlated
#
# Usage: source() this file in R. Adjust parameters (seed, n_assoc, alpha) as desired.

set.seed(123)  # for reproducibility; change if you want different random choices

# ---- Parameters ----
marker_names <- c("CD8", "CD4", "PD-L1", "FOXP3", "vimentin",
                  "Ki67", "panCK", "CD68", "CD31", "GranzymeB")

# Strength of correlation: larger alpha -> stronger positive effect of survival on intensity
alpha <- 0.8

# How many (marker, cell_type) associations to create
n_assoc <- 10

# Coordinate bounds for spatial simulation
x_min <- 0; x_max <- 100
y_min <- 0; y_max <- 100

# ---- Sample-level summary including survival times and counts ----
set.seed(6868)
sample_summary <- data.frame(
  SampleID = c("OS001","OS002","OS003","OS004","OS005","OS006","OS007","OS008","OS009","OS010"),
  Survival_time = rnbinom(10, size = 5, mu = 5),  # years
  Tumor_cells = rnbinom(10, size = 10, mu = 200),
  CD8_Tcells = rnbinom(10, size = 10, mu = 90),
  CD4_Tcells = rnbinom(10, size = 10, mu = 100),
  Macrophage = rnbinom(10, size = 10, mu = 75),
  Stromal_cells = rnbinom(10, size = 10, mu = 100),
  Fibroblast = rnbinom(10, size = 10, mu = 80),
  Endothelial = rnbinom(10, size = 10, mu = 35),
  NK_cell = rnbinom(10, size = 10, mu = 12),
  Dendritic_cell = rnbinom(10, size = 5, mu = 7),
  stringsAsFactors = FALSE
)

# Identify cell-type columns
cell_type_cols <- setdiff(names(sample_summary), c("SampleID", "Survival_time"))

# ---- Prepare baseline parameters for simulation ----
# Baseline intensity per marker (vector)
baseline_marker <- setNames(rgamma(length(marker_names), shape = 2, scale = 0.5),  
                            # positive, skewed,
                            marker_names)

# Small cell-type offset to make marker means vary by cell-type
celltype_offset_vals <- runif(length(cell_type_cols), min = -0.15, max = 0.15)
names(celltype_offset_vals) <- cell_type_cols

# Scale survival_time to 0..1 for modulation
surv <- sample_summary$Survival_time
surv_scaled <- (surv - min(surv)) / (max(surv) - min(surv))

# ---- Choose random (marker, cell_type) pairs that will correlate with survival ----
all_pairs <- expand.grid(marker = marker_names, cell_type = cell_type_cols,
                         stringsAsFactors = FALSE)
assoc_idx <- sample(seq_len(nrow(all_pairs)), n_assoc)
assoc_pairs <- all_pairs[assoc_idx, ]
assoc_pairs$alpha <- alpha  # same strength for all; could be varied per pair

# Save associations so user can inspect which pairs were used
write.csv(assoc_pairs, "survival_associations.csv", row.names = FALSE)

message("Selected survival-associated (marker, cell_type) pairs:")
print(assoc_pairs)

# ---- Simulation loop: create a file per sample ----
for (i in seq_len(nrow(sample_summary))) {
  sample_id <- sample_summary$SampleID[i]
  s_scaled <- surv_scaled[i]
  
  # Prepare empty list to collect rows (faster than repeated rbind)
  rows <- vector("list", sum(as.numeric(sample_summary[i, cell_type_cols])))
  row_counter <- 0
  cell_num <- 1
  
  for (cell_type in cell_type_cols) {
    n_cells <- sample_summary[[cell_type]][i]
    if (n_cells <= 0) next
    
    for (j in seq_len(n_cells)) {
      row_counter <- row_counter + 1
      # Base mean for markers for this cell (marker baseline + celltype offset)
      # We'll then add the survival effect for associated pairs and noise.
      # Use reproducible small jitter per cell-type/marker if desired.
      Xcoord <- runif(1, x_min, x_max)
      Ycoord <- runif(1, y_min, y_max)
      
      # Prepare marker values
      marker_vals <- numeric(length(marker_names))
      names(marker_vals) <- marker_names
      
      for (m in marker_names) {
        base_mean <- baseline_marker[m] + celltype_offset_vals[cell_type]
        
        # Check if this (marker,cell_type) is an associated pair
        is_assoc <- any(assoc_pairs$marker == m & assoc_pairs$cell_type == cell_type)
        if (is_assoc) {
          beta <- assoc_pairs$alpha[assoc_pairs$marker == m & assoc_pairs$cell_type == cell_type]
          mean_int <- base_mean + beta * s_scaled
        } else {
          mean_int <- base_mean
        }
        # Add per-cell noise
        val <- rnorm(1, mean = mean_int, sd = 0.12)
        # Clip to plausible range (0 to 2)
        val <- min(max(val, 0), 2)
        marker_vals[m] <- round(val, 3)
      }
      
      # Build the row
      row <- c(
        CellID = paste0("C", sprintf("%05d", cell_num)),
        Type = cell_type,
        X = round(Xcoord, 2),
        Y = round(Ycoord, 2),
        as.list(marker_vals)
      )
      
      rows[[row_counter]] <- row
      cell_num <- cell_num + 1
    } # end cells loop
  } # end cell_type loop
  
  # Convert list of rows to data.frame
  df <- do.call(rbind.data.frame, rows)
  # rbind.data.frame produced factors for strings on some setups; coerce properly:
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  # Convert numeric columns back to numeric
  df$X <- as.numeric(df$X); df$Y <- as.numeric(df$Y)
  for (m in marker_names) df[[m]] <- as.numeric(df[[m]])
  
  # Write to CSV
  outname <- paste0(sample_id, "_cells.csv")
  write.csv(df, outname, row.names = FALSE)
  message("Wrote ", outname, " (", nrow(df), " cells).")
}

message("Simulation complete. See survival_associations.csv for the (marker, cell_type) pairs correlated with survival.")


# 
# set.seed(1234)
# rnbinom(10, size = 20, mu = 200)