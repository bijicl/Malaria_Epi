# Load required packages
library(dplyr)
library(ggplot2)
library(forcats)
library(broom)
# Load libraries
pacman::p_load(readr, dplyr, segmented, ggplot2)
pacman::p_load(dplyr, zoo, readr)


# --- STEP 1: Read and prepare data ---
# Replace with your actual file path if using CSV
malaria_data <- read.csv("aggregated_data_malaria_clinical.csv")
malaria_data <- malaria_data %>%
  mutate(
    year_month_clean   = as.yearmon(as.character(year_month), format = "%Y-%m"),
    year_month_numeric = as.numeric(year_month_clean)
  )
#  Step 2: Filter valid rows
malaria_clean <- malaria_data %>%
  filter(!is.na(year_month_numeric) & !is.na(MalariaResults_P))

# Step 3: Fit base linear model
lm_model <- lm(MalariaResults_P ~ year_month_numeric, data = malaria_clean)

#  Step 4: Pick a psi guess inside the year range
range(malaria_clean$year_month_numeric)
# (for example: between 2010.5 and 2020)

# 🔀 Step 5: Fit joinpoint regression
seg_model <- segmented(lm_model, seg.Z = ~year_month_numeric,
                       psi = list(year_month_numeric = c(2015.5)))  # example guess

#  Step 6: View summary
summary(seg_model)

#  Step 7: Plot the fitted trend and joinpoint
ggplot(malaria_clean, aes(x = year_month_clean, y = MalariaResults_P)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_line(aes(y = fitted(seg_model)), color = "blue", linewidth = 1.2) +
  geom_vline(xintercept = as.yearmon(seg_model$psi[, "Est."]), linetype = "dashed", color = "firebrick") +
  annotate("text", x = as.yearmon(seg_model$psi[, "Est."]),
           y = max(malaria_clean$MalariaResults_P),
           label = paste("Joinpoint @", format(as.yearmon(seg_model$psi[, "Est."]), "%Y-%m")),
           vjust = 1) +
  theme_minimal() +
  labs(title = "Joinpoint Regression: Malaria Positivity Over Time",
       x = "Year-Month", y = "Malaria Positive Cases")


