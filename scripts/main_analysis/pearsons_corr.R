m <- TMS_metadata %>% group_by(Age_group, mouse.id)
m <- m %>% filter(Senescence.group == "Senescence.pos")
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(Age_group, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a)
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*10000
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
write.csv(a, "senescence_norm_to_10000.csv")

a <- read.csv("senescence_norm_to_10000.csv")
a <- a[,-1]
c <- read.csv("CD4_subsets_without_zero_counts.csv")
c <- c[,-1]
c <- read.csv("CD4_subsets_with_zero_counts.csv")
c <- c[,-1]
a$No_of_senescence_cells <- a$No_of_cells
a$No_of_cells <- NULL
a$Normalize_senescence_cells_per_10000 <- a$normalize
a$normalize <- NULL
d <- merge(a, c)
cor(d$Normalize_senescence_cells_per_10000, d$No_of_cells, method = "pearson", use = "complete.obs")
ggscatter(d, x = "Normalize_senescence_cells_per_10000", y = "percent", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "senescence per 10000", ylab = "cd4 subset number")
d <- d %>% filter(subsets == "Naive")
d <- d %>% filter(subsets == "Naive_Isg15")
d <- d %>% filter(subsets == "rTregs")
d <- d %>% filter(subsets == "aTregs")
d <- d %>% filter(subsets == "Cytotoxic")
d <- d %>% filter(subsets == "TEM")
d <- d %>% filter(subsets == "Exhausted")
d <- d %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
d <- d %>% filter(Age_group.x == "Young")
d <- d %>% filter(Age_group == "Old")

#CD4_subsets_without_zero_counts
m <- CD4_T_metadata %>% group_by(Age_group, mouse.id, subsets)
m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
a <- merge(m,a)
a$percent <- "0"
for (n in 1:length(a[,5])) {
  a$percent[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$percent <- as.numeric(as.vector(a$percent))
write.csv(a, "CD4_subsets_without_zero_counts.csv")


#CD4_subsets_with_zero_counts
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a)
a$percent <- "0"
for (n in 1:length(a[,5])) {
  a$percent[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$percent <- as.numeric(as.vector(a$percent))
write.csv(a, "CD4_subsets_with_zero_counts.csv")


combat <- d[,]
d <- round(cor(d), 2)


#CD4/CD8 correlation with senescence
m <- CD8_T_metadata %>% group_by(Age_group, mouse.id)
#m <- m %>% summarise()
m <- m %>% summarise(No_of_CD8_T_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)
a <- a %>% filter(subsets != "Not T")
a <- a %>% summarise(No_of_CD4_T_cells = n())
a <- merge(m,a)
a$CD4_to_CD8 <- "0"
for (n in 1:length(a[,5])) {
  a$CD4_to_CD8[n] <- a$No_of_CD4_T_cells[n]/a$No_of_CD8_T_cells[n]
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$CD4_to_CD8 <- as.numeric(as.vector(a$CD4_to_CD8))

c <- read.csv("senescence_norm_to_10000.csv")
c <- c[,-1]
c <- merge(c,a)
c$Normalize_senescence_cells_per_10000 <- c$normalize
c$normalize <- NULL
c$Normalize_senescence_cells_per_10000 <- as.numeric(as.vector(c$Normalize_senescence_cells_per_10000))
cor(c$Normalize_senescence_cells_per_10000, c$CD4_to_CD8, method = "pearson", use = "complete.obs")
ggscatter(c, x = "Normalize_senescence_cells_per_10000", y = "CD4_to_CD8", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "senescence per 10000", ylab = "CD4/CD8")
c <- c %>% filter(Age_group == "Old")
c <- c %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )

#corr age with sen
TMS_metadata <- read.csv("senescence_per_total_cells.csv")
m <- TMS_metadata %>% group_by(Age_group, age, mouse.id)
m <- m %>% filter(Senescence.group == "Senescence.pos")
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(Age_group, age, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a)
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*10000
}
a$age <- plyr::mapvalues(x = a$age, from = c("1m", "3m", "18m", "21m", "24m", "30m"),
                                            to = c("1", "3", "18", "21", "24", "30"))
a$age <- as.numeric(as.vector(a$age))
a$normalize <- as.numeric(as.vector(a$normalize))
cor(a$age, a$normalize, method = "pearson", use = "complete.obs")
ggscatter(a, x = "age", y = "normalize", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "age", ylab = "senescence per 10000", panel.labs = c("Old", "Young"))
library(corrplot)
ggqqplot(a$age, ylab = "age")
ggqqplot(a$normalize, ylab = "normalize")
rcorr(a$age, a$normalize, type = c("pearson","spearman"))
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}
cor_3 <- rcorr(as.matrix(a))

my_cor_matrix <- flat_cor_mat(cor_3$r, cor_3$P)