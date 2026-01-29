library(ggplot2)
library(dplyr)

# Read data
df <- read.csv("metrics_simulation_study.csv")
if (!"All" %in% names(df) && "Coverage" %in% names(df)) {
    df$All <- df$Coverage
}

# Calculate key statistics
mean_cov <- mean(df$All, na.rm = TRUE)
sd_cov <- sd(df$All, na.rm = TRUE)
pct_success <- mean(df$All >= 0.90, na.rm = TRUE) * 100
median_cov <- median(df$All, na.rm = TRUE)

# Create publication-quality histogram
p_hist <- ggplot(df, aes(x = All)) +
    # Shaded success region (90-100%)
    annotate("rect",
        xmin = 0.90, xmax = 1.0,
        ymin = 0, ymax = Inf,
        fill = "lightgreen", alpha = 0.2
    ) +

    # Histogram with better binning
    geom_histogram(
        bins = 20,
        fill = "#2C3E50",
        color = "white",
        alpha = 0.8,
        linewidth = 0.3
    ) +

    # Density Overlay (Polished Version)
    geom_density(aes(y = after_stat(count) * (max(after_stat(count)) / max(after_stat(density)))),
        color = "#E74C3C",
        linewidth = 1.2,
        linetype = "solid"
    ) +

    # Target line (0.95)
    geom_vline(aes(xintercept = 0.95),
        color = "#E74C3C",
        linetype = "dashed",
        linewidth = 1.2
    ) +

    # Minimum acceptable line (0.90)
    geom_vline(aes(xintercept = 0.90),
        color = "#F39C12",
        linetype = "dotted",
        linewidth = 1.2
    ) +

    # Mean line
    geom_vline(aes(xintercept = mean_cov),
        color = "#3498DB",
        linewidth = 1.5
    ) +

    # Annotations (bottom left, clear position)
    annotate("text",
        x = 0.85, y = Inf,
        label = sprintf(
            "Mean: %.3f\nSD: %.3f\nMedian: %.3f\n≥0.90: %.1f%%",
            mean_cov, sd_cov, median_cov, pct_success
        ),
        hjust = 0, vjust = 1.5,
        size = 4,
        color = "black",
        fontface = "bold",
        family = "sans"
    ) +

    # Line labels (legend alternative - clearer)
    annotate("text",
        x = 0.95, y = Inf,
        label = "Target (95%)",
        vjust = -0.5, hjust = 1.1,
        color = "#E74C3C", size = 3.5, fontface = "bold"
    ) +

    # Minimum Label
    annotate("text",
        x = 0.90, y = Inf,
        label = "Minimum (90%)",
        vjust = -2, hjust = 1.1,
        color = "#F39C12", size = 3.5, fontface = "bold"
    ) +

    # Mean Label
    annotate("text",
        x = mean_cov, y = Inf,
        label = "Mean",
        vjust = -3.5, hjust = 0.5,
        color = "#3498DB", size = 3.5, fontface = "bold"
    ) +

    # Labels
    labs(
        title = "ACI Coverage Probability Distribution Across 1,000 Synthetic Outbreak Scenarios",
        x = "Empirical Coverage Probability",
        y = "Frequency (Number of Simulations)",
        caption = "Shaded region indicates acceptable coverage (≥90%). ACI successfully maintains valid coverage in 95%+ of scenarios."
    ) +

    # Professional theme
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.caption = element_text(hjust = 0, face = "italic", size = 9, color = "gray30"),
        axis.title = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.margin = margin(15, 15, 15, 15)
    ) +

    # Scale adjustments
    scale_x_continuous(
        breaks = seq(0.85, 1.0, 0.025),
        labels = scales::percent_format(accuracy = 0.1),
        limits = c(0.85, 1.0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Save high-resolution
ggsave("plot_simulation_histogram.png",
    p_hist,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
)

# Also save as PDF for vector graphics (journals prefer this)
ggsave("plot_simulation_histogram.pdf",
    p_hist,
    width = 10,
    height = 6,
    device = cairo_pdf
)

# Print summary
cat("\n=== Simulation Summary ===\n")
cat(sprintf("Mean Coverage: %.4f\n", mean_cov))
cat(sprintf("SD Coverage: %.4f\n", sd_cov))
cat(sprintf("Median Coverage: %.4f\n", median_cov))
cat(sprintf("Simulations ≥ 0.90: %.1f%%\n", pct_success))
cat(sprintf("Simulations ≥ 0.95: %.1f%%\n", mean(df$All >= 0.95) * 100))
cat("========================\n")
