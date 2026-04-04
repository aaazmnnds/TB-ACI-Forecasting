library(tidyverse)
library(lubridate)

results_dir <- "results"
input_file <- file.path(results_dir, "metrics_ACI_rolling_eval.csv")
output_file <- file.path(results_dir, "poro_covid_metrics.csv")

df <- read.csv(input_file) %>% mutate(Date = as.Date(DateObj))

# COVID period: March 2020 - December 2021
covid_data <- df %>%
    filter(
        Date >= as.Date("2020-03-01"),
        Date <= as.Date("2021-12-31")
    )

poro_metrics <- data.frame(
    Coverage = mean(covid_data$Covered),
    Mean_Width = mean(covid_data$Width)
)

write.csv(poro_metrics, output_file, row.names = FALSE)
cat("Poro COVID metrics saved to:", output_file, "\n")
