set.seed(43)

data_wo_CAF_v2  <- matrix(c(
  "H", 0.043, 0.0012, 0.0072, 0.23, 0.000022, 0.000000054,
  "H", 0.038, 0.00098, 0.0072, 0.23, 0.000022, 0.000000054,
  "H", 0.05, 0.0012, 0.0072, 0.23, 0.000022, 0.000000054,
  "H", 0.04, 0.001, 0.0072, 0.23, 0.000022, 0.000000054,
  "H", 0.035, 0.0012, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.034, 0.0023, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.03, 0.00098, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.035, 0.001, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.031, 0.0019, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.025, 0.00048, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.073, 0.00057, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.058, 0.0008, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.04, 0.0011, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.069, 0.0014, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.061, 0.0011, 0.0072, 0.23, 0.000022, 0.000000054
), nrow = 15, ncol = 7, byrow = TRUE)

df <- as.data.frame(data_wo_CAF_v2, stringsAsFactors = FALSE)

colnames(df) <- c("E2", 
                  "k1_hat","alpha1","beta","db","ER0","E2ER0")

numeric_cols <- c("k1_hat","alpha1","beta","db","ER0","E2ER0")
df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)

# Convert treatment to factor
df$E2 <- factor(df$E2, levels = c("H", "M", "L"),
                       labels = c("H", "M", "L"))

library(tidyr)

df_long <- pivot_longer(
  df,
  cols = -E2,
  names_to = "parameter",
  values_to = "value"
)

library(ggplot2)

ggplot(df_long, aes(x = E2, y = value, fill = E2)) +
  geom_boxplot() +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(title = "Parameter Distributions by E2 supply",
       x = "E2",
       y = "Value") +
  theme_bw()

### k1hat and alpha1 #################
data_wo_CAF_v2  <- matrix(c(
  "H", 0.043, 0.0012, 0.0072, 0.23, 0.000022, 0.000000054,
  "H", 0.038, 0.00098, 0.0072, 0.23, 0.000022, 0.000000054,
  "H", 0.05, 0.0012, 0.0072, 0.23, 0.000022, 0.000000054,
  "H", 0.04, 0.001, 0.0072, 0.23, 0.000022, 0.000000054,
  "H", 0.035, 0.0012, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.034, 0.0023, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.03, 0.00098, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.035, 0.001, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.031, 0.0019, 0.0072, 0.23, 0.000022, 0.000000054,
  "M", 0.025, 0.00048, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.073, 0.00057, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.058, 0.0008, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.04, 0.0011, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.069, 0.0014, 0.0072, 0.23, 0.000022, 0.000000054,
  "L", 0.061, 0.0011, 0.0072, 0.23, 0.000022, 0.000000054
), nrow = 15, ncol = 7, byrow = TRUE)

df <- as.data.frame(data_wo_CAF_v2[,1:3], stringsAsFactors = FALSE)

colnames(df) <- c("E2", 
                  "k1hat","alpha1")

numeric_cols <- c("k1hat","alpha1")
df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)

# Convert treatment to factor
df$E2 <- factor(df$E2, levels = c("H", "M", "L"),
                labels = c("H", "M", "L"))

library(tidyr)

df_long <- pivot_longer(
  df,
  cols = -E2,
  names_to = "parameter",
  values_to = "value"
)

library(ggplot2)

ggplot(df_long, aes(x = E2, y = value, fill = E2)) +
  geom_boxplot() +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(title = "Parameter Distributions by E2 supply",
       x = "E2",
       y = "Value") +
  theme_bw()
################
data_with_CAF <- matrix(c(
  "H", 0.17, 13.04, 90.25, 0.000049,
  "H", 0.13, 13.01, 28.85, 0.00005,
  "H", 0.23, 12.94, 137.14, 0.000048,
  "H", 0.13, 13.42, 106.55, 0.00005,
  "H", 0.2, 13.45, 281.28, 0.000048,
  "M", 0.082, 13.79, 15.77, 0.000052,
  "M", 0.11, 13.67, 32.42, 200.09,
  "M", 0.094, 13.49, 35.89, 0.000051,
  "M", 0.21, 12.39, 18.33, 158.77,
  "M", 0.081, 13.63, 16.62, 0.000051,
  "L", 0.18, 13.12, 73.88, 0.000049,
  "L", 0.17, 12.74, 17.93, 57.85,
  "L", 0.22, 12.49, 54.11, 0.000049,
  "L", 0.27, 11.67, 93.94, 0.000045,
  "L", 0.2, 12.44, 43.65, 0.000049
), nrow = 15, ncol = 5, byrow = TRUE)

df2 <- as.data.frame(data_with_CAF, stringsAsFactors = FALSE)

colnames(df2) <- c("E2", 
                  "k2", "alpha2", "d3", "alpha3")

numeric_cols <- c("k2", "alpha2", "d3", "alpha3")
df2[numeric_cols] <- lapply(df2[numeric_cols], as.numeric)

# Convert treatment to factor
df2$E2 <- factor(df2$E2, levels = c("H", "M", "L"),
                labels = c("H", "M", "L"))

library(tidyr)

df2_long <- pivot_longer(
  df2,
  cols = -E2,
  names_to = "parameter",
  values_to = "value"
)

library(ggplot2)

ggplot(df2_long, aes(x = E2, y = value, fill = E2)) +
  geom_boxplot() +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(title = "Parameter Distributions by E2 supply",
       x = "E2",
       y = "Value") +
  theme_bw()
