library(ggplot2)
df_opt <- read.csv("../../results/optimization_history.csv")
p_opt <- ggplot(df_opt, aes(x = trial_number, y = value)) + 
  geom_line(color="lightgray") + geom_point(aes(color = value), size=3) +
  geom_point(data = df_opt[which.min(df_opt$value),], color="red", size=5, shape=1) +
  scale_color_viridis_c(alpha=0.8) +
  labs(title="Bayesian Hyperparameter Optimization History (50 Trials)",
       x="Trial Number", y="Objective Value L(lambda,W) = 100|PICP - 0.95| + MPIW", color="Obj") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("../../Submission_Package_JBI/Fig2_Optuna_History.png", p_opt, width=8, height=5, dpi=300, bg="white")
ggsave("../../Fig2_Optuna_History.png", p_opt, width=8, height=5, dpi=300, bg="white")
