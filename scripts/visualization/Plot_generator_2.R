library(ggpubr)
library(ggplot2)

CD8_T_metadata <- read.csv("CD8_T_metedata_TMS.csv")
CD4_T_metadata <- read.csv("CD4_T_metedata_TMS.csv")
TMS_metadata <- read.csv("metedata_TMS.csv")
b <- read.csv("CD4_T_subsets.csv")
b <- b[,-1]

#Number of CD8
m <- CD8_T_metadata %>% group_by(Age_group_2, mouse.id)
#m <- m %>% summarise()
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(Age_group_2, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*10000
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 1000), 
                                                                             breaks = c(0, 200, 400, 600, 800, 1000)) +
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(980)) + 
  theme_classic() + scale_shape_manual(values=1:23)

write.csv(a, "CD8_T_per_total_cells.csv")

#Number of CD4
m <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id)
#m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(Age_group_2, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*10000
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 600), 
                                                                             breaks = c(0, 200, 400, 600)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(550)) + 
  theme_classic()+ scale_shape_manual(values=1:23)

write.csv(a, "CD4_T_per_total_cells.csv")

#Number of Senescence
m <- TMS_metadata %>% group_by(Age_group, mouse.id)
m <- m %>% filter(Senescence.group == "Senescence.pos")
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(Age_group, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a[is.na(a)] <- 0
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*10000
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group, color = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 200), breaks = c(0, 50, 100, 150, 200)) +
  stat_compare_means(aes(label = paste0("p = ", ..p.adj..)), label.y = 150, label.x = 1.25) + 
  theme_classic()+ scale_shape_manual(values=1:22) + xlab("Age group") + 
  ylab("Normalize senescence cells per 10000 cells per mouse") +  
  ggtitle("Senescence cells amount per total cells vs age group") +
  scale_color_manual(values=c("#64b4b6", "#897457")) + 
  scale_size_manual(values=c(3,3)) + 
  theme(text = element_text(size = 20))

write.csv(a, "senescence_per_total_cells.csv")

a <- a %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 200), breaks = c(0, 200)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(100, 150, 200)) + 
  theme_classic()+ scale_shape_manual(values=1:22)

write.csv(a, "senescence_per_total_cells_without_fucked_mice.csv")

#CD4/CD8
m <- CD8_T_metadata %>% group_by(Age_group_2, mouse.id)
#m <- m %>% summarise()
m <- m %>% summarise(No_of_CD8_T_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id)
#a <- a %>% filter(subsets != "Not T")
a <- a %>% summarise(No_of_CD4_T_cells = n())
a <- merge(m,a, all = TRUE)
a$CD4_to_CD8 <- "0"
for (n in 1:length(a[,5])) {
  a$CD4_to_CD8[n] <- a$No_of_CD4_T_cells[n]/a$No_of_CD8_T_cells[n]
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$CD4_to_CD8 <- as.numeric(as.vector(a$CD4_to_CD8))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = CD4_to_CD8)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = CD4_to_CD8), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 2), breaks = c(0, 0.5, 1, 1.5, 2)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(2)) + 
  theme_classic()+ scale_shape_manual(values=1:23)+ 
  ylab("CD4/CD8 per mouse") + 
  xlab("Age group 2") + 
  ggtitle("CD4/CD8 per mouse vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035")) 

write.csv(a, "CD4_to_CD8.csv")

a <- a %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
ggplot(a, aes(x = Age_group, y = CD4_to_CD8)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = CD4_to_CD8), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 2), breaks = c(0, 0.5, 1, 1.5, 2)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(2)) + 
  theme_classic()+ scale_shape_manual(values=1:22)

write.csv(a, "CD4_to_CD8_without_fucked_mice.csv")


#CD4 subsets Naive without zero counts
m <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id, subsets)
#m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id) %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE) %>% filter(subsets == "Naive")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 90), breaks = c(0, 90)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(90)) + 
  theme_classic()+ scale_shape_manual(values=1:22)+ 
  ylab("Percent of Naive subset out of CD4 T cells per mouse") + 
  xlab("Age group 2") + 
  ggtitle("Percent of Naive subset out of CD4 T cells per mouse vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035")) 
write.csv(a, "naive_per_age_without_zero_counts.csv")

#CD4 subsets Naive with zero counts
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Naive")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 90), breaks = c(0, 90)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(90)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "naive_per_age.csv")

#CD4 subsets Naive without mice with total CD4 cells under 40: 3-M-5/6, 3-M-7/8, 3-M-9, 18-M-52, 18-M-53, 24-M-61, 
#30-M-2
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Naive")
a <- a %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 90), breaks = c(0, 90)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(90)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "naive_per_age_without_mice_with_total_CD4_cells_under_40.csv")


#CD4 subsets Naive_Isg15 without zero counts
m <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id, subsets)
#m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id) %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE) %>% filter(subsets == "Naive_Isg15")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 35), breaks = c(0, 35)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(30)) + 
  theme_classic()+ scale_shape_manual(values=1:22)+ 
  ylab("Percent of Naive subset out of CD4 T cells per mouse") + 
  xlab("Age group 2") + 
  ggtitle("Percent of Naive subset out of CD4 T cells per mouse vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035")) 
write.csv(a, "naive_isg15_per_age_without_zero_counts.csv")

#CD4 subsets Naive_Isg15 with zero counts
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Naive_Isg15")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 5), breaks = c(0, 5)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(5)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "naive_isg15_per_age.csv")

#CD4 subsets Naive_Isg15 without mice with total CD4 cells under 40: 3-M-5/6, 3-M-7/8, 3-M-9, 18-M-52, 18-M-53, 24-M-61, 
#30-M-2
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Naive_Isg15")
a <- a %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 5), breaks = c(0, 5)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(5)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "naive_isg15_per_age_without_mice_with_total_CD4_cells_under_40.csv")

#CD4 subsets rTregs without zero counts
m <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id, subsets)
#m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id) %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE) %>% filter(subsets == "rTregs")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 20), breaks = c(0, 20)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(19)) + 
  theme_classic()+ scale_shape_manual(values=1:22)+ 
  ylab("Percent of rTregs subset out of CD4 T cells per mouse") + 
  xlab("Age group 2") + 
  ggtitle("Percent of rTregs subset out of CD4 T cells per mouse vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035")) 
write.csv(a, "rTregs_per_age_without_zero_counts.csv")


#CD4 subsets rTregs with zero counts
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "rTregs")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 35), breaks = c(0, 3500)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(35)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "rTregs_per_age.csv")


#CD4 subsets rTregs without mice with total CD4 cells under 40: 3-M-5/6, 3-M-7/8, 3-M-9, 18-M-52, 18-M-53, 24-M-61, 
#30-M-2
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "rTregs")
a <- a %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 35), breaks = c(0, 3500)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(35)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "rTregs_per_age_without_mice_with_total_CD4_cells_under_40.csv")


#CD4 subsets aTregs without zero counts
m <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id, subsets)
#m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id) %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE) %>% filter(subsets == "aTregs")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 50), breaks = c(0, 50)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(50)) + 
  theme_classic()+ scale_shape_manual(values=1:22)+ 
  ylab("Percent of aTregs subset out of CD4 T cells per mouse") + 
  xlab("Age group 2") + 
  ggtitle("Percent of aTregs subset out of CD4 T cells per mouse vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035")) 
write.csv(a, "aTregs_per_age_without_zero_counts.csv")


#CD4 subsets aTregs with zero counts
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "aTregs")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 30), breaks = c(0, 30)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(20, 25, 30)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "aTregs_per_age.csv")


#CD4 subsets aTregs without mice with total CD4 cells under 40: 3-M-5/6, 3-M-7/8, 3-M-9, 18-M-52, 18-M-53, 24-M-61, 
#30-M-2
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "aTregs")
a <- a %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 30), breaks = c(0, 30)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(20, 25, 30)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "aTregs_per_age_without_mice_with_total_CD4_cells_under_40.csv")


#CD4 subsets Cytotoxic without zero counts
m <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id, subsets)%>% filter(subsets == "Cytotoxic")
#m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id) %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)# %>% filter(subsets == "Cytotoxic")
for (n in 1:length(a[,5])) {
  if (is.na(a$No_of_cells[n])) a$No_of_cells[n] <- 0
}
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 70), breaks = c(0, 70)) +
  stat_compare_means( method = "t.test", comparisons = my_comparisons, label.y = c(70)) + 
  theme_classic()+ scale_shape_manual(values=1:22)+ 
  ylab("Percent of Cytotoxic subset out of CD4 T cells per mouse") + 
  xlab("Age group 2") + 
  ggtitle("Percent of Cytotoxix subset out of CD4 T cells per mouse vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035"))
write.csv(a, "cytotoxix_per_age_without_zero_counts.csv")


#CD4 subsets Cytotoxic with zero counts
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Cytotoxic")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 50), breaks = c(0, 50)) +
  stat_compare_means( method = "t.test", comparisons = my_comparisons, label.y = c(40, 45, 50)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "cytotoxix_per_age.csv")

#CD4 subsets Cytotoxic without mice with total CD4 cells under 40: 3-M-5/6, 3-M-7/8, 3-M-9, 18-M-52, 18-M-53, 24-M-61, 
#30-M-2
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Cytotoxic")
a <- a %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 50), breaks = c(0, 50)) +
  stat_compare_means( method = "t.test", comparisons = my_comparisons, label.y = c(40, 45, 50)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "cytotoxix_per_age_without_mice_with_total_CD4_cells_under_40.csv")


#CD4 subsets TEM without zero counts
m <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id, subsets)
m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "TEM")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 50), breaks = c(0, 50)) +
  stat_compare_means( method = "t.test", comparisons = my_comparisons, label.y = c(50)) + 
  theme_classic()+ scale_shape_manual(values=1:22)+ 
  ylab("Percent of TEM subset out of CD4 T cells per mouse") + 
  xlab("Age group 2") + 
  ggtitle("Percent of TEM subset out of CD4 T cells per mouse vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035"))
write.csv(a, "tem_per_age_without_zero_counts.csv")


#CD4 subsets TEM with zero counts
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "TEM")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 50), breaks = c(0, 50)) +
  stat_compare_means( method = "t.test", comparisons = my_comparisons, label.y = c(40, 45, 50)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "tem_per_age.csv")


#CD4 subsets TEM without mice with total CD4 cells under 40: 3-M-5/6, 3-M-7/8, 3-M-9, 18-M-52, 18-M-53, 24-M-61, 
#30-M-2
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "TEM")
a <- a %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 50), breaks = c(0, 50)) +
  stat_compare_means( method = "t.test", comparisons = my_comparisons, label.y = c(40, 45, 50)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "tem_per_age_without_mice_with_total_CD4_cells_under_40.csv")


#CD4 subsets Exhausted without zero counts
m <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id, subsets)
m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Exhausted")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 40), breaks = c(0, 40)) +
  stat_compare_means( method = "t.test", comparisons = my_comparisons, label.y = c(30, 35, 40)) + 
  theme_classic()+ scale_shape_manual(values=1:22)+ 
  ylab("Percent of Exhausted subset out of CD4 T cells per mouse") + 
  xlab("Age group 2") + 
  ggtitle("Percent of Exhausted subset out of CD4 T cells per mouse vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035"))
write.csv(a, "exhausted_per_age_without_zero_counts.csv")


#CD4 subsets Exhausted with zero counts
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Exhausted")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 40), breaks = c(0, 40)) +
  stat_compare_means( method = "t.test", comparisons = my_comparisons, label.y = c(30, 35, 40)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "exhausted_per_age.csv")


#CD4 subsets Exhausted without mice with total CD4 cells under 40: 3-M-5/6, 3-M-7/8, 3-M-9, 18-M-52, 18-M-53, 24-M-61, 
#30-M-2
m <- b
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% 
  summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Exhausted")
a <- a %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" )
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$No_of_total_CD4_T_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 40), breaks = c(0, 40)) +
  stat_compare_means( method = "t.test", comparisons = my_comparisons, label.y = c(30, 35, 40)) + 
  theme_classic()+ scale_shape_manual(values=1:22)
write.csv(a, "exhausted_per_age_without_mice_with_total_CD4_cells_under_40.csv")

#CD4 subsets Naive
m <- CD4_T_metadata %>% group_by(Age_group, mouse.id, subsets)
m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% summarise(total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Naive")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_CD4_T_cells[n]*100
}
a$subsets <- factor(a$subsets, levels = c("Naive", "Naive_Isg15", "rTregs", "aTregs", "TEM", "Cytotoxic", "Exhausted"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Naive", "Naive_Isg15"), c("Naive", "rTregs"), c("Naive", "aTregs"), c("Naive", "TEM"), 
                        c("Naive", "Cytotoxic"), c("Naive", "Exhausted"), c("Naive_Isg15", "rTregs"), 
                        c("Naive_Isg15", "aTregs"), c("Naive_Isg15", "TEM"), c("Naive_Isg15", "Cytotoxic"),
                        c("Naive_Isg15", "Exhausted"), c("rTregs", "aTregs"), c("rTregs", "TEM"), 
                        c("rTregs", "Cytotoxic"), c("rTregs", "Exhausted"), c("aTregs", "TEM"), 
                        c("Naive", "Cytotoxic"), c("Naive", "Exhausted"),)
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 90), breaks = c(0, 90)) +
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(80, 85, 90)) + 
  theme_classic()
write.csv(a, "naive_per_age.csv")

#CD4 subsets Naive
m <- CD4_T_metadata %>% group_by(Age_group, mouse.id, subsets)
m <- m %>% filter(subsets != "Not T")
m <- m %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group, mouse.id)%>% filter(subsets != "Not T") %>% summarise(total_CD4_T_cells = n())
#normalized total cells number per mice
a <- merge(m,a) %>% filter(subsets == "Naive")
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_CD4_T_cells[n]*100
}
a$subsets <- factor(a$subsets, levels = c("Naive", "Naive_Isg15", "rTregs", "aTregs", "TEM", "Cytotoxic", "Exhausted"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Naive", "Naive_Isg15"), c("Naive", "rTregs"), c("Naive", "aTregs"), c("Naive", "TEM"), 
                        c("Naive", "Cytotoxic"), c("Naive", "Exhausted"), c("Naive_Isg15", "rTregs"), 
                        c("Naive_Isg15", "aTregs"), c("Naive_Isg15", "TEM"), c("Naive_Isg15", "Cytotoxic"),
                        c("Naive_Isg15", "Exhausted"), c("rTregs", "aTregs"), c("rTregs", "TEM"), 
                        c("rTregs", "Cytotoxic"), c("rTregs", "Exhausted"), c("aTregs", "TEM"), 
                        c("Naive", "Cytotoxic"), c("Naive", "Exhausted"),)
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 90), breaks = c(0, 90)) +
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(80, 85, 90)) + 
  theme_classic()
write.csv(a, "naive_per_age.csv")