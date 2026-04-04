library(tidyverse)
library(lubridate)

# --------------------------------------------------------------------------------
# Brazil SRAG Data Preparation (Simplified for External Validation)
# Focusing on consistent municipality-level historical data (2002-2023)
# --------------------------------------------------------------------------------

# 1. Configuration
set.seed(2026)
# Robust Path Resolution
find_results_dir <- function() {
    if (dir.exists("results")) {
        return("results")
    }
    if (dir.exists("../../results")) {
        return("../../results")
    }
    if (dir.exists("../results")) {
        return("../results")
    }
    return("results") # Fallback
}
results_dir <- find_results_dir()
cat(sprintf("Using results directory: %s\n", normalizePath(results_dir)))

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Robust Input Resolution
find_input_file <- function(filename) {
    if (file.exists(filename)) {
        return(filename)
    }
    if (file.exists(file.path("../../", filename))) {
        return(file.path("../../", filename))
    }
    if (file.exists(file.path("../", filename))) {
        return(file.path("../", filename))
    }
    if (file.exists(file.path(results_dir, filename))) {
        return(file.path(results_dir, filename))
    }
    return(filename)
}

# Historical data (only consistent source)
historical_csv <- find_input_file("brazil_sivep_gripe.csv")
output_csv <- file.path(results_dir, "brazil_sivep_gripe_extended_2025.csv") # Kept same name for script compatibility

cat("Loading historical data (2002-2023)...\n")
df <- read.csv(historical_csv) %>%
    mutate(Date = as.Date(Date))

# Standardize columns
if ("Population_est" %in% names(df)) {
    df <- df %>% rename(Population = Population_est)
}

# Ensure final formatting (Standardizing to Poro Schema)
df_final <- df %>%
    arrange(Date) %>%
    mutate(
        Year = lubridate::year(Date),
        Month = lubridate::month(Date, label = TRUE, abbr = FALSE),
        Incidence_per_100k = (Count / Population) * 100000
    )

cat(sprintf("Writing cleaned dataset to %s...\n", output_csv))
write.csv(df_final, output_csv, row.names = FALSE)

# Visualize series and structural break
library(ggplot2)

df_final$Date <- as.Date(df_final$Date)

p <- ggplot(df_final, aes(x = Date, y = Count)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_vline(
        xintercept = as.Date("2020-03-01"),
        linetype = "dashed", color = "red", linewidth = 1
    ) +
    annotate("text",
        x = as.Date("2020-03-01"), y = max(df_final$Count) * 0.85,
        label = "COVID-19 Onset (March 2020)", color = "red", angle = 90, vjust = -0.5, fontface = "bold"
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))) +
    labs(
        title = "Brazil SRAG Monthly Notifications (2002-2023)",
        subtitle = "Consolidated municipal-level historical surveillance data",
        x = "Year", y = "Monthly Case Count"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14)
    )

ggsave(file.path(results_dir, "plot_brazil_structural_break.png"), p, width = 10, height = 5, bg = "white", dpi = 300)
ggsave(file.path(results_dir, "plot_brazil_structural_break.pdf"), p, width = 10, height = 5)
cat(sprintf("\nStructural break visualization saved to %s\n", results_dir))
cat("Done!\n")
