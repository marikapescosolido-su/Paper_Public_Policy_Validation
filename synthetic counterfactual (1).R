##############################
# 1. Load Libraries
##############################
library(dplyr)
library(lubridate)
library(Synth)
library(ggplot2)
library(zoo)
library(tidyr)
library(fixest)
library(did)
library(sdid)

##############################
# 2. Load Raw Panel
##############################
df_full <- read.csv("C:/Users/Jinho Kim/Desktop/GATech/2025 Fall/PUBP 3042/pm25_panel.csv",
               stringsAsFactors = FALSE)

# Convert to proper Date
df_full$date <- ymd(df_full$date)

##############################
# 3. Collapse duplicates
# One row per city per month
##############################
# Collapse duplicates to monthly level in the FULL data
df_full <- df_full %>%
  mutate(year_month = format(date, "%Y-%m")) %>%
  group_by(country, city, year_month) %>%
  summarise(pm25 = mean(pm25, na.rm = TRUE), .groups = "drop") %>%
  mutate(date = as.Date(paste0(year_month, "-01")))

##############################
# 4. Restrict to balanced window
# Ensures perfect balance for SCM
##############################
# ---------- Subset for SCM (balanced pre + around treatment) ----------
df_sc <- df_full %>%
  filter(date >= as.Date("2019-06-01"),
         date <= as.Date("2021-01-01")) %>%
  arrange(city, date)

df_sc$year_month <- format(df_sc$date, "%Y-%m")
df_sc$unit_id <- as.numeric(as.factor(df_sc$city))
months_vec_sc <- sort(unique(df_sc$year_month))
df_sc$time_id <- match(df_sc$year_month, months_vec_sc)

##############################
# 6. Define treatment info
##############################
treated_city <- "Karachi"
treated_id <- df_sc$unit_id[df_sc$city == treated_city][1]
donor_ids   <- unique(df_sc$unit_id[df_sc$city != treated_city])

treat_month <- "2021-01"
treat_time  <- which(months_vec_sc == treat_month)
pre_period  <- which(months_vec_sc < treat_month)
full_period_sc <- 1:length(months_vec_sc)

##############################
# 7. Synth Data Prep
##############################
df_sc_df <- as.data.frame(df_sc)

dataprep.out <- dataprep(
  foo = df_sc_df,
  predictors = "pm25",
  dependent  = "pm25",
  unit.variable = "unit_id",
  time.variable = "time_id",
  treatment.identifier  = treated_id,
  controls.identifier   = donor_ids,
  time.predictors.prior = pre_period,
  time.optimize.ssr     = pre_period,
  time.plot             = full_period_sc
)

##############################
# 8. Run SYNTH
##############################
synth.out <- synth(dataprep.out)

############################################
# 9. Compute Synthetic Karachi for ALL months
############################################
# Get weights and donor cities
W <- synth.out$solution.w              # matrix of weights
control_ids <- as.numeric(rownames(W)) # these are unit_ids of donors in df_sc
W_vec <- as.numeric(W)

# Map unit_ids back to city names
id_to_city <- df_sc %>%
  distinct(unit_id, city)

donor_cities <- id_to_city$city[match(control_ids, id_to_city$unit_id)]

# Build wide panel from full data
panel_wide <- df_full %>%
  filter(city %in% c(treated_city, donor_cities)) %>%
  arrange(date) %>%
  distinct(city, year_month, .keep_all = TRUE) %>%
  select(city, year_month, date, pm25) %>%
  tidyr::pivot_wider(
    id_cols = c(year_month, date),
    names_from = city,
    values_from = pm25
  ) %>%
  arrange(date)

# Restrict to months where all needed donors & Karachi are present
needed_cols <- c(treated_city, donor_cities)
panel_wide <- panel_wide %>%
  filter(if_all(all_of(needed_cols), ~ !is.na(.)))

# Compute synthetic Karachi for all these months
Y1_full <- panel_wide[[treated_city]]
X_donors <- as.matrix(panel_wide[, donor_cities])  # columns order matches donor_cities
Y0_full <- as.numeric(X_donors %*% W_vec)

results <- data.frame(
  date      = panel_wide$date,
  year_month = panel_wide$year_month,
  actual    = Y1_full,
  synthetic = Y0_full
)

results$period <- ifelse(results$date < as.Date("2021-01-01"), "Pre", "Post")

##############################
# 10. Plot Actual vs Synthetic (Raw) for full time
##############################
ggplot(results, aes(x = date)) +
  geom_line(aes(y = actual,    colour = "Karachi (Actual)"), linewidth = 1.1) +
  geom_line(aes(y = synthetic, colour = "Synthetic Karachi"), linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-01"), linetype = "dotted") +
  scale_color_manual(values = c("Karachi (Actual)" = "black",
                                "Synthetic Karachi" = "red")) +
  labs(
    title = "Synthetic Control: Karachi PM2.5 (Full Period)",
    x = "Date",
    y = "PM2.5 (µg/m³)",
    colour = ""
  ) +
  theme_minimal(base_size = 14)

##############################
# 11. STL (Seasonal Adjustment)
##############################
ts_actual    <- ts(results$actual,    frequency = 12)
ts_synthetic <- ts(results$synthetic, frequency = 12)

stl_actual    <- stl(ts_actual, s.window="periodic")
stl_synthetic <- stl(ts_synthetic, s.window="periodic")

# Seasonally adjusted (remove seasonal component)
results$actual_stl    <- as.numeric(ts_actual    - stl_actual$time.series[,"seasonal"])
results$synthetic_stl <- as.numeric(ts_synthetic - stl_synthetic$time.series[,"seasonal"])

##############################
# 12. Plot STL-adjusted Comparison
##############################

ggplot(results, aes(x = date)) +
  geom_line(aes(y = actual_stl,    colour = "Karachi STL-adjusted"), linewidth = 1.1) +
  geom_line(aes(y = synthetic_stl, colour = "Synthetic STL-adjusted"),
            linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-01"), linetype = "dotted") +
  scale_color_manual(values = c("Karachi STL-adjusted" = "blue",
                                "Synthetic STL-adjusted" = "orange")) +
  labs(
    title = "STL-Adjusted Synthetic Control: Karachi PM2.5",
    x = "Date",
    y = "PM2.5 (µg/m³, seasonally adjusted)",
    colour = ""
  ) +
  theme_minimal(base_size = 14)
##############################
# 13. Gap Plots (Causal Effect)
##############################
# Raw gap
results$gap <- results$actual - results$synthetic

ggplot(results, aes(x = date, y = gap)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = as.Date("2021-01-01"), linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Euro-V Effect on Karachi PM2.5 (Actual - Synthetic)",
    x = "Date",
    y = "PM2.5 difference (µg/m³)"
  ) +
  theme_minimal(base_size = 14)

# STL adjusted gap
results$gap_stl <- results$actual_stl - results$synthetic_stl

ggplot(results, aes(x = date, y = gap_stl)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = as.Date("2021-01-01"), linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Euro-V Effect on Karachi PM2.5 (STL-Adjusted)",
    x = "Date",
    y = "Seasonally-adjusted gap (Actual - Synthetic)"
  ) +
  theme_minimal(base_size = 14)

##############################
# 14. Perform DiD
##############################
# Create DiD variables
# First: reshape into long format (actual vs synthetic as two "units")
did_df <- results %>%
  select(date, year_month, actual, synthetic) %>%
  pivot_longer(cols = c(actual, synthetic),
               names_to = "unit_type",
               values_to = "pm25") %>%
  mutate(
    treated = ifelse(unit_type == "actual", 1, 0),
    post    = ifelse(date >= as.Date("2021-01-01"), 1, 0),
    time_id = as.integer(as.factor(year_month)),
    unit_id = ifelse(unit_type == "actual", 1, 0)
  )
##############################
# A. Two-Way Fixed Effects DiD
##############################

did_twfe <- feols(
  pm25 ~ treated * post | unit_id + time_id,
  data = did_df
)

summary(did_twfe)
##############################

