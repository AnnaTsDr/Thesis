library(tidyverse)

dat <- as.data.frame(c)
dat$mouse.id <- NULL
dat$Age_group_2 <- NULL

dat_old <- c %>% filter(Age_group_2 == "Old")
dat_old <- as.data.frame(dat_old)
dat_old$mouse.id <- NULL
dat_old$Age_group_2 <- NULL

round(cor(dat), digits = 2)

library(corrplot)
corrplot(cor(dat), method = "number", type = "upper")


cor(c$Cytotoxic, c$Senescence, method = "pearson")

corrplot2 <- function(data, method = "pearson", sig.level = 0.05, order = "original", diag = FALSE, type = "upper",
                      tl.srt = 90, number.font = 1, number.cex = 1, mar = c(0, 0, 0, 0)) {
  data_incomplete <- data
  data <- data[complete.cases(data), ]
  mat <- cor(data, method = method)
  cor.mtest <- function(mat, method) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = method)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  p.mat <- cor.mtest(data, method = method)
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(mat,
           method = "color", col = col(200), number.font = number.font,
           mar = mar, number.cex = number.cex,
           type = type, order = order,
           addCoef.col = "black", # add correlation coefficient
           tl.col = "black", tl.srt = tl.srt, # rotation of text labels
           # combine with significance level
           p.mat = p.mat, sig.level = sig.level, insig = "blank",
           # hide correlation coefficients on the diagonal
           diag = diag
  )
}

# edit from here
corrplot2(
  data = dat,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75
)

library("PerformanceAnalytics")
chart.Correlation(dat, histogram=TRUE, pch=19)

chart.Correlation(dat_old, histogram=TRUE, pch=19)

res <- cor(dat)
col<- colorRampPalette(c("#2700D1", "white", "#CC0C00"))(20)
heatmap(x = res, col = col, symm = TRUE)

res2 <- cor(dat_old)
heatmap(x = res2, col = col, symm = TRUE)
