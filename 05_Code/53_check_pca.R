suppressMessages(library(readr))
ex <- as.matrix(read.csv("results/00_Data_Prep/corrected_expr.csv",
                         row.names = 1, check.names = FALSE))
pc <- prcomp(t(ex), scale. = TRUE)
imp <- summary(pc)$importance[2, 1:2]
cat("PC1:", round(100 * imp[1], 1), "%  PC2:", round(100 * imp[2], 1), "%\n")
cat("first 5 sample coordinates (PC1, PC2):\n")
print(round(pc$x[1:5, 1:2], 3))
