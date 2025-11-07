# Step 0: Load and install required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, broom, forcats, forestplot, grid)

# Step 1: Input malaria data
malaria <- data.frame(
  Age_Group = c("<1 year", "1-4 years", "5-8 years", "9-14 years", ">14 years"),
  PF = c(212, 583, 564, 640, 4172),
  PV = c(58, 30, 29, 71, 670),
  Mixed = c(5, 3, 5, 1, 24)
)

# Step 2: Calculate totals and non-infection counts
malaria <- malaria %>%
  mutate(
    Total = Mixed + PF + PV,
    Non_PF = Total - PF,
    Non_PV = Total - PV,
    Non_Mixed = Total - Mixed
  )

# Step 3: Function to expand and model
expand_and_model <- function(outcome, non_outcome, label) {
  # Always use dplyr::select to avoid conflicts
  infected <- malaria %>%
    dplyr::select(Age_Group, infection_count = all_of(outcome)) %>%
    tidyr::uncount(weights = infection_count) %>%
    mutate(y = 1)
  non_infected <- malaria %>%
    dplyr::select(Age_Group, non_count = all_of(non_outcome)) %>%
    tidyr::uncount(weights = non_count) %>%
    mutate(y = 0)
  long_df <- bind_rows(infected, non_infected)
  
  long_df$Age_Group <- factor(long_df$Age_Group,
                              levels = c(">14 years", "9-14 years", "5-8 years", "1-4 years", "<1 year"))
  
  model <- glm(y ~ Age_Group, data = long_df, family = binomial)
  results <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(Infection = label,
           Age_Group = gsub("Age_Group", "", term),
           Age_Group = gsub("C\\(|\\)\\[T\\.|\\]", "", Age_Group))
  return(results)
}

# Step 4: Run models for PF, PV, Mixed
pf_results <- expand_and_model("PF", "Non_PF", "PF")
pv_results <- expand_and_model("PV", "Non_PV", "PV")
mixed_results <- expand_and_model("Mixed", "Non_Mixed", "Mixed")
# Step 5: Combine and view
all_results <- bind_rows(pf_results, pv_results, mixed_results) %>%
  dplyr::select(Infection, Age_Group, estimate, conf.low, conf.high, p.value) %>%
  arrange(Infection, Age_Group)

# Optional: print results
print(all_results)

# Reorder age groups chronologically
age_order <- c("<1 year", "1-4 years", "5-8 years", "9-14 years", ">14 years")

# Reorder Infection levels explicitly
infection_order <- c("PF", "PV", "Mixed")

# Modify final output with custom order
all_results <- all_results %>%
  mutate(
    Age_Group = factor(Age_Group, levels = age_order),
    Infection = factor(Infection, levels = infection_order)
  ) %>%
  arrange(Infection, Age_Group)

# Optional: view results
print(all_results)

# Create color vectors based on significance
box_colors <- ifelse(all_results$p.value < 0.05, "red", "gray")  # Adjust colors as desired
line_colors <- ifelse(all_results$p.value < 0.05, "darkred", "gray40")

# Define font style
fpFont <- grid::gpar(fontfamily = "Times", fontsize = 10)

# Label matrix for display (excluding Significance column)
labeltext <- cbind(
  c("Age Group", as.character(all_results$Age_Group)),
  c("Species Type", as.character(all_results$Infection)),
  c("OR [95% CI]", sprintf("%.2f [%.2f–%.2f]", all_results$estimate, all_results$conf.low, all_results$conf.high)),
  c("P-value", sprintf("%.3f", all_results$p.value))
)

# Forestplot
forestplot::forestplot(
  labeltext = labeltext,
  mean  = c(NA, all_results$estimate),
  lower = c(NA, all_results$conf.low),
  upper = c(NA, all_results$conf.high),
  is.summary = c(TRUE, rep(FALSE, nrow(all_results))),
  xlab = "Odds Ratio with 95% Confidence Interval",
  title = "Forest Plot of Malaria Odds by Age Group and Species Type",
  xticks = exp(seq(log(0.1), log(10), length.out = 6)),
  zero = 1,
  grid = TRUE,
  boxsize = 0.2,
  col = forestplot::fpColors(box = box_colors,
                             line = line_colors,
                             summary = "black"),
  txt_gp = forestplot::fpTxtGp(
    label = list(fpFont, fpFont),
    ticks = fpFont,
    xlab = fpFont,
    title = fpFont
  ),
  graph.pos = 4,
  align = "llrr",
  clip = c(0.1, 10)
)

# Optional: Legend for color coding
legend("topright", 
       legend = c("Significant (p < 0.05)", "Non-significant"), 
       col = c("red", "gray"), 
       pch = 15, 
       bty = "n", 
       cex = 0.8)

# Trial 3 Mutation Type vs. Non-Type Logistic Regression
# Load required libraries
if (!require("forestplot")) install.packages("forestplot")
library(forestplot)
library(grid)
library(dplyr)
library(tidyr)
library(broom)

# Step 1: Input mutation data
mutation_data <- data.frame(
  Type = c("Point Mutation", "Deletion", "Duplication", "Frameshift", "Stop-Gain", "Other"),
  Count = c(21900, 763, 340, 0, 311, 2579)
)

# Step 2: Calculate total and non-type counts
mutation_data <- mutation_data %>%
  mutate(
    Total = sum(Count),
    Non_Count = Total - Count
  )

# Step 3: Function to expand and model each mutation type
expand_and_model <- function(type_label) {
  case_count <- mutation_data %>% filter(Type == type_label) %>% pull(Count)
  non_count <- mutation_data %>% filter(Type == type_label) %>% pull(Non_Count)
  
  # Skip if case_count is zero
  if (case_count == 0) return(NULL)
  
  # Expand case rows
  case_df <- data.frame(
    Type = type_label,
    y = 1
  ) %>% slice(rep(1, case_count))
  
  # Expand non-case rows
  non_df <- data.frame(
    Type = paste0("Not_", type_label),
    y = 0
  ) %>% slice(rep(1, non_count))
  
  # Combine
  long_df <- bind_rows(case_df, non_df)
  
  # Fit logistic regression
  model <- glm(y ~ Type, data = long_df, family = binomial)
  
  # Extract results
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(Mutation = type_label)
}

# Step 4: Run models for all mutation types with non-zero counts
valid_types <- mutation_data %>% filter(Count > 0) %>% pull(Type)
results_list <- lapply(valid_types, expand_and_model)
all_results <- bind_rows(results_list) %>%
  select(Mutation, estimate, conf.low, conf.high, p.value) %>%
  mutate(
    OR_CI = sprintf("%.2f [%.2f–%.2f]", estimate, conf.low, conf.high),
    Pval = sprintf("%.3f", p.value)
  )

# Step 5: Prepare forest plot labels
labeltext <- rbind(
  c("Mutation Type", "OR [95% CI]", "P-value"),
  cbind(
    as.character(all_results$Mutation),
    all_results$OR_CI,
    all_results$Pval
  )
)

# Step 6: Forest plot styling
box_colors <- ifelse(all_results$p.value < 0.05, "red", "gray80")
line_colors <- ifelse(all_results$p.value < 0.05, "darkred", "gray40")
fpFont <- gpar(fontfamily = "Times", fontsize = 10)

# Step 7: Forest plot
forestplot(
  labeltext = labeltext,
  mean  = c(NA, all_results$estimate),
  lower = c(NA, all_results$conf.low),
  upper = c(NA, all_results$conf.high),
  is.summary = c(TRUE, rep(FALSE, nrow(all_results))),
  xlab = "Odds Ratio (Mutation Type vs. All Others)",
  title = "Forest Plot: Mutation Type Significance",
  xticks = exp(seq(log(0.1), log(10), length.out = 6)),
  zero = 1,
  grid = TRUE,
  boxsize = 0.2,
  col = fpColors(box = box_colors,
                 line = line_colors,
                 summary = "black"),
  txt_gp = fpTxtGp(
    label = list(fpFont, fpFont),
    ticks = fpFont,
    xlab = fpFont,
    title = fpFont
  ),
  graph.pos = 3,
  align = "llr",
  clip = c(0.1, 10)
)

# Optional: Legend
legend("topright",
       legend = c("Significant (p < 0.05)", "Non-significant"),
       col = c("red", "gray80"),
       pch = 15,
       bty = "n",
       cex = 0.8)