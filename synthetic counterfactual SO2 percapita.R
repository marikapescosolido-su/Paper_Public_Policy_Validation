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
# library(sdid)

##############################
# 2. Load Raw SO2 Panel (wide) and reshape to long
##############################
library(dplyr)
library(tidyr)
library(lubridate)

df_full <- read.csv(
  "C:/Users/daga1/OneDrive/Documents/PUBP3042 Project/EDGAR_SO2_m_2000_2022_v2/SO2_monthly_percapita_wide.csv",
  stringsAsFactors = FALSE
)
# Check columns
str(df_full)

# Month columns in the wide file
month_cols <- c("Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec")
month_cols_pc <- paste0(month_cols, "_kg_per_person")

# 1) Wide -> Long (per-capita SO2)
df_full <- df_full %>%
  tidyr::pivot_longer(
    cols = all_of(month_cols_pc),
    names_to = "MonthCol",
    values_to = "so2"        # still call it 'so2' so the rest of the code doesn't change
  ) %>%
  mutate(
    MonthName = sub("_kg_per_person", "", MonthCol)  # "Jan_kg_per_person" -> "Jan"
  )

# 2) Map month names to numbers
month_map <- setNames(1:12, month.abb)  # Jan -> 1, Feb -> 2, ...

df_full <- df_full %>%
  mutate(
    month = month_map[MonthName],
    # create Date from Year + month
    date  = as.Date(paste(year, month, "01", sep = "-")),
    year_month = format(date, "%Y-%m"),
    name = country
  )

# sanity check
head(df_full)

# Convert to proper Date
df_full$date <- ymd(df_full$date)

##############################
# 3. Collapse duplicates
# One row per country per month
##############################
df_full <- df_full %>%
  mutate(year_month = format(date, "%Y-%m")) %>%
  group_by(country, year_month) %>%
  summarise(so2 = mean(so2, na.rm = TRUE), .groups = "drop") %>%
  mutate(date = as.Date(paste0(year_month, "-01")))

##############################
# 4. Restrict to balanced window
##############################
df_sc <- df_full %>%
  filter(date >= as.Date("2008-01-01"),
         date <= as.Date("2021-01-01")) %>%
  arrange(country, date)

df_sc$year_month <- format(df_sc$date, "%Y-%m")

# Unit IDs are now countries instead of cities
df_sc$unit_id <- as.numeric(as.factor(df_sc$country))

months_vec_sc <- sort(unique(df_sc$year_month))
df_sc$time_id <- match(df_sc$year_month, months_vec_sc)

##############################
# 6. Define treatment info
##############################
treated_country <- "Pakistan"
treated_id <- df_sc$unit_id[df_sc$country == treated_country][1]
donor_ids   <- unique(df_sc$unit_id[df_sc$country != treated_country])

treat_month <- "2021-01"
treat_time  <- which(months_vec_sc == treat_month)
pre_period  <- which(months_vec_sc < treat_month)
full_period_sc <- 1:length(months_vec_sc)

##############################
# 7. Synth Data Prep (NOx per capita)
##############################
df_sc_df <- as.data.frame(df_sc)

dataprep.out <- dataprep(
  foo = df_sc_df,
  predictors = "so2",
  dependent  = "so2",
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
# 9. Compute Synthetic Pakistan for ALL months
############################################
W <- synth.out$solution.w
control_ids <- as.numeric(rownames(W))
W_vec <- as.numeric(W)

id_to_country <- df_sc %>%
  distinct(unit_id, country)

donor_countries <- id_to_country$country[match(control_ids, id_to_country$unit_id)]

panel_wide <- df_full %>%
  filter(country %in% c(treated_country, donor_countries)) %>%
  arrange(date) %>%
  distinct(country, year_month, .keep_all = TRUE) %>%
  select(country, year_month, date, so2) %>%
  tidyr::pivot_wider(
    id_cols = c(year_month, date),
    names_from = country,
    values_from = so2
  ) %>%
  arrange(date)

needed_cols <- c(treated_country, donor_countries)
panel_wide <- panel_wide %>%
  filter(if_all(all_of(needed_cols), ~ !is.na(.)))

Y1_full <- panel_wide[[treated_country]]
X_donors <- as.matrix(panel_wide[, donor_countries])
Y0_full <- as.numeric(X_donors %*% W_vec)

results <- data.frame(
  date       = panel_wide$date,
  year_month = panel_wide$year_month,
  actual     = Y1_full,
  synthetic  = Y0_full
)

results$period <- ifelse(results$date < as.Date("2021-01-01"), "Pre", "Post")

##############################
# 10. Plot Actual vs Synthetic (Raw) for full time
##############################
ggplot(results, aes(x = date)) +
  geom_line(aes(y = actual,    colour = "Pakistan (Actual)"), linewidth = 1.1) +
  geom_line(aes(y = synthetic, colour = "Synthetic Pakistan"), linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-01"), linetype = "dotted") +
  scale_color_manual(values = c("Pakistan (Actual)" = "black",
                                "Synthetic Pakistan" = "red")) +
  labs(
    title = "Synthetic Control: Pakistan SO2 per capita (Full Period)",
    x = "Date",
    y = "SO2 per capita (kg/person/month)",
    colour = ""
  ) +
  theme_minimal(base_size = 14)

##############################
# 11. STL (Seasonal Adjustment)
##############################
ts_actual    <- ts(results$actual,    frequency = 12)
ts_synthetic <- ts(results$synthetic, frequency = 12)

stl_actual    <- stl(ts_actual, s.window="periodic", robust = TRUE)
stl_synthetic <- stl(ts_synthetic, s.window="periodic", robust = TRUE)

results$actual_stl    <- as.numeric(ts_actual    - stl_actual$time.series[,"seasonal"])
results$synthetic_stl <- as.numeric(ts_synthetic - stl_synthetic$time.series[,"seasonal"])

##############################
# 12. Plot STL-adjusted Comparison
##############################
ggplot(results, aes(x = date)) +
  geom_line(aes(y = actual_stl,    colour = "Pakistan STL-adjusted"), linewidth = 1.1) +
  geom_line(aes(y = synthetic_stl, colour = "Synthetic STL-adjusted"),
            linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-01"), linetype = "dotted") +
  scale_color_manual(values = c("Pakistan STL-adjusted" = "blue",
                                "Synthetic STL-adjusted" = "orange")) +
  labs(
    title = "STL-Adjusted Synthetic Control: Pakistan SO2 per capita",
    x = "Date",
    y = "SO2 per capita (kg/person/month, seasonally adjusted)",
    colour = ""
  ) +
  theme_minimal(base_size = 14)

##############################
# 13. Gap Plots (Causal Effect)
##############################
results$gap <- results$actual - results$synthetic

ggplot(results, aes(x = date, y = gap)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = as.Date("2021-01-01"), linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Euro-V Effect on Pakistan SO2 per capita (Actual - Synthetic)",
    x = "Date",
    y = "SO2 per capita difference (kg/person/month)"
  ) +
  theme_minimal(base_size = 14)

results$gap_stl <- results$actual_stl - results$synthetic_stl

ggplot(results, aes(x = date, y = gap_stl)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = as.Date("2021-01-01"), linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Euro-V Effect on Pakistan SO2 per capita (STL-Adjusted)",
    x = "Date",
    y = "Seasonally-adjusted gap (kg/person/month)"
  ) +
  theme_minimal(base_size = 14)

##############################
# 14. Perform DiD (unchanged, just interprets units as per capita)
##############################
did_df <- results %>%
  select(date, year_month, actual, synthetic) %>%
  pivot_longer(cols = c(actual, synthetic),
               names_to = "unit_type",
               values_to = "so2") %>%
  mutate(
    treated = ifelse(unit_type == "actual", 1, 0),
    post    = ifelse(date >= as.Date("2021-01-01"), 1, 0),
    time_id = as.integer(as.factor(year_month)),
    unit_id = ifelse(unit_type == "actual", 1, 0)
  )

did_twfe <- feols(
  so2 ~ treated * post | unit_id + time_id,
  data = did_df
)
summary(did_twfe)

##############################
# 15. Placebo test: in-time placebo before Jan 2021
##############################
placebo_date <- as.Date("2020-01-01")

results_pre <- results %>%
  filter(date < as.Date("2021-01-01"))

did_df_placebo <- results_pre %>%
  select(date, year_month, actual, synthetic) %>%
  pivot_longer(cols = c(actual, synthetic),
               names_to = "unit_type",
               values_to = "so2") %>%
  mutate(
    treated = ifelse(unit_type == "actual", 1, 0),
    post    = ifelse(date >= placebo_date, 1, 0),
    time_id = as.integer(as.factor(year_month)),
    unit_id = ifelse(unit_type == "actual", 1, 0)
  )

did_twfe_placebo <- feols(
  so2 ~ treated * post | unit_id + time_id,
  data = did_df_placebo
)
summary(did_twfe_placebo)

ggplot(results_pre, aes(x = date)) +
  geom_line(aes(y = actual,    colour = "Pakistan (Actual)"), linewidth = 1.1) +
  geom_line(aes(y = synthetic, colour = "Synthetic Pakistan"), linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = placebo_date, linetype = "dotted") +
  scale_color_manual(values = c("Pakistan (Actual)" = "black",
                                "Synthetic Pakistan" = "red")) +
  labs(
    title = "Placebo Test: Raw Actual vs Synthetic (Pre-treatment Only)",
    subtitle = paste("Fake treatment at", placebo_date),
    x = "Date",
    y = "SO2 per capita (kg/person/month)",
    colour = ""
  ) +
  theme_minimal(base_size = 14)

ggplot(results_pre, aes(x = date)) +
  geom_line(aes(y = actual_stl,    colour = "Pakistan STL-adjusted"), linewidth = 1.1) +
  geom_line(aes(y = synthetic_stl, colour = "Synthetic STL-adjusted"),
            linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = placebo_date, linetype = "dotted") +
  scale_color_manual(values = c("Pakistan STL-adjusted" = "blue",
                                "Synthetic STL-adjusted" = "orange")) +
  labs(
    title = "Placebo Test: STL-Adjusted Actual vs Synthetic (Pre-treatment Only)",
    subtitle = paste("Fake treatment at", placebo_date),
    x = "Date",
    y = "SO2 per capita (kg/person/month, seasonally adjusted)",
    colour = ""
  ) +
  theme_minimal(base_size = 14)

results_pre$gap_stl <- results_pre$actual_stl - results_pre$synthetic_stl

did_df_placebo_stl <- results_pre %>%
  select(date, year_month, actual_stl, synthetic_stl) %>%
  pivot_longer(cols = c(actual_stl, synthetic_stl),
               names_to = "unit_type",
               values_to = "so2") %>%
  mutate(
    treated = ifelse(unit_type == "actual_stl", 1, 0),
    post    = ifelse(date >= placebo_date, 1, 0),
    time_id = as.integer(as.factor(year_month)),
    unit_id = ifelse(unit_type == "actual_stl", 1, 0)
  )

did_twfe_placebo_stl <- feols(
  so2 ~ treated * post | unit_id + time_id,
  data = did_df_placebo_stl
)
summary(did_twfe_placebo_stl)