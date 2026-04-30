#Number of Senescence 3 age groups
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
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old", "Supercentenarian"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"), c("Supercentenarian", "Old"), c("Young", "Supercentenarian"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group, color = Age_group)) + 
  stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 250), breaks = c(0, 50, 100, 150, 200, 250)) +
  stat_compare_means(aes(label = paste0("p = ", ..p.adj..), method = "t.test"), comparisons = my_comparisons, 
                     label.y = c(200, 220, 240), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group") + 
  ylab("Normalize senescence cells per 10000 cells per mouse") +  
  ggtitle("Senescence cells amount per total cells vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035", "#331900")) + 
  scale_size_manual(values=c(3,3,3)) + 
  theme(text = element_text(size = 20))

#write.csv(a, "senescence_per_total_cells.csv")

#Number of Senescence 2 age groups
m <- TMS_metadata %>% group_by(Age_group_2, mouse.id)
m <- m %>% filter(Senescence.group == "Senescence.pos")
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(Age_group_2, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a[is.na(a)] <- 0
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*10000
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"))
ggplot(a, aes(x = Age_group_2, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + 
  stat_summary(aes(x = Age_group_2, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 250), breaks = c(0, 50, 100, 150, 200, 250)) +
  stat_compare_means(aes(label = paste0("p = ", ..p.adj..), method = "t.test"), comparisons = my_comparisons, 
                     label.y = c(240), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group") + 
  ylab("Normalize senescence cells per 10000 cells per mouse") +  
  ggtitle("Senescence cells amount per total cells vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035")) + 
  scale_size_manual(values=c(3,3)) + 
  theme(text = element_text(size = 20))

#Number of Senescence age
m <- TMS_metadata %>% group_by(age, mouse.id)
m <- m %>% filter(Senescence.group == "Senescence.pos")
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(age, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a[is.na(a)] <- 0
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*10000
}
a$age <- factor(a$age, levels = c("1m", "3m", "18m", "21m", "24m", "30m"))
a$normalize <- as.numeric(as.vector(a$normalize))
ggplot(a, aes(x = age, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = age, color = age)) + 
  stat_summary(aes(x = age, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 250), breaks = c(0, 50, 100, 150, 200, 250)) +
  stat_compare_means(aes(label = paste0("ANOVA p = ", ..p.adj..), method = "anova"), 
                     label.y = c(240), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group") + 
  ylab("Normalize senescence cells per 10000 cells per mouse") +  
  ggtitle("Senescence cells amount per total cells vs age group") +
  scale_color_manual(values=c("#00A9CC", "#65EFFF", "#F2DACD", "#CC9B7A", "#996035", "#331900")) + 
  scale_size_manual(values=c(3,3,3,3,3,3)) + 
  theme(text = element_text(size = 20))

#Number of p16 age
m <- TMS_metadata %>% group_by(age, mouse.id)
m <- m %>% filter(p16.group == "p16.pos")
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(age, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a[is.na(a)] <- 0
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*10000
}
a$age <- factor(a$age, levels = c("1m", "3m", "18m", "21m", "24m", "30m"))
a$normalize <- as.numeric(as.vector(a$normalize))
ggplot(a, aes(x = age, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = age, color = age)) + 
  stat_summary(aes(x = age, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 250), breaks = c(0, 50, 100, 150, 200, 250)) +
  stat_compare_means(aes(label = paste0("ANOVA p = ", ..p.adj..), method = "anova"), 
                     label.y = c(240), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group") + 
  ylab("Normalize p16 positive cells per 10000 cells per mouse") +  
  ggtitle("p16 positive cells cells amount per total cells vs age group") +
  scale_color_manual(values=c("#00A9CC", "#65EFFF", "#F2DACD", "#CC9B7A", "#996035", "#331900")) + 
  scale_size_manual(values=c(3,3,3,3,3,3)) + 
  theme(text = element_text(size = 20))

#Number of p21 age
m <- TMS_metadata %>% group_by(age, mouse.id)
m <- m %>% filter(p21.group == "p21.pos")
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(age, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a[is.na(a)] <- 0
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*10000
}
a$age <- factor(a$age, levels = c("1m", "3m", "18m", "21m", "24m", "30m"))
a$normalize <- as.numeric(as.vector(a$normalize))
ggplot(a, aes(x = age, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = age, color = age)) + 
  stat_summary(aes(x = age, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 10000), breaks = c(0, 2000, 4000, 6000, 8000, 10000)) +
  stat_compare_means(aes(label = paste0("ANOVA p = ", ..p.adj..), method = "anova"), 
                     label.y = c(9000), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group") + 
  ylab("Normalize p21 positive cells per 10000 cells per mouse") +  
  ggtitle("p21 positive cells amount per total cells vs age group") +
  scale_color_manual(values=c("#00A9CC", "#65EFFF", "#F2DACD", "#CC9B7A", "#996035", "#331900")) + 
  scale_size_manual(values=c(3,3,3,3,3,3)) + 
  theme(text = element_text(size = 20))

#Number of cells per age
m <- TMS_metadata %>% group_by(age, mouse.id)
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(age) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a[is.na(a)] <- 0
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*100
}
a$age <- factor(a$age, levels = c("1m", "3m", "18m", "21m", "24m", "30m"))
a$normalize <- as.numeric(as.vector(a$normalize))
a$No_of_cells <- as.numeric(as.vector(a$No_of_cells))
ggplot(a, aes(x = age, y = No_of_cells)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = age, color = age)) + 
  stat_summary(aes(x = age, y = No_of_cells), fun = median,  fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 30000), breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000)) +
  stat_compare_means(aes(label = paste0("ANOVA p = ", ..p.adj..), method = "anova"), 
                     label.y = c(28000), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group") + 
  ylab("Normalize p21 positive cells per 10000 cells per mouse") +  
  ggtitle("p21 positive cells amount per total cells vs age group") +
  scale_color_manual(values=c("#00A9CC", "#65EFFF", "#F2DACD", "#CC9B7A", "#996035", "#331900")) + 
  scale_size_manual(values=c(3,3,3,3,3,3)) + 
  theme(text = element_text(size = 20))

#Number of cells per Age_group
m <- TMS_metadata %>% group_by(Age_group, mouse.id)
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(Age_group) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a[is.na(a)] <- 0
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*100
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old", "Supercentenarian"))
a$normalize <- as.numeric(as.vector(a$normalize))
a$No_of_cells <- as.numeric(as.vector(a$No_of_cells))
ggplot(a, aes(x = Age_group, y = No_of_cells)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group, color = Age_group)) + 
  stat_summary(aes(x = Age_group, y = No_of_cells), fun = median,  fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 30000), breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000)) +
  stat_compare_means(aes(label = paste0("ANOVA p = ", ..p.adj..), method = "anova"), 
                     label.y = c(28000), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group") + 
  ylab("Normalize p21 positive cells per 10000 cells per mouse") +  
  ggtitle("p21 positive cells amount per total cells vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035", "#331900")) + 
  scale_size_manual(values=c(3,3,3,3,3,3)) + 
  theme(text = element_text(size = 20))

#Number of cells per Age_group_2
m <- TMS_metadata %>% group_by(Age_group_2, mouse.id)
m <- m %>% summarise(No_of_cells = n())
a <- TMS_metadata %>% group_by(Age_group_2) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a[is.na(a)] <- 0
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_cells[n]/a$total_cells[n]*100
}
a$Age_group_2 <- factor(a$Age_group_2, levels = c("Young", "Old", "Supercentenarian"))
a$normalize <- as.numeric(as.vector(a$normalize))
a$No_of_cells <- as.numeric(as.vector(a$No_of_cells))
ggplot(a, aes(x = Age_group_2, y = No_of_cells)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + 
  stat_summary(aes(x = Age_group_2, y = No_of_cells), fun = median,  fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 30000), breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000)) +
  stat_compare_means(aes(label = paste0("ANOVA p = ", ..p.adj..), method = "anova"), 
                     label.y = c(28000), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group") + 
  ylab("Normalize p21 positive cells per 10000 cells per mouse") +  
  ggtitle("p21 positive cells amount per total cells vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035", "#331900")) + 
  scale_size_manual(values=c(3,3,3,3,3,3)) + 
  theme(text = element_text(size = 20))

#Number of cells per mouse per Age_group_2
m <- TMS_metadata %>% group_by(Age_group_2, mouse.id)
m <- m %>% summarise(No_of_cells = n())

m$Age_group_2 <- factor(m$Age_group_2, levels = c("Young", "Old"))
m$No_of_cells <- as.numeric(as.vector(m$No_of_cells))

ggplot(m, aes(x = Age_group_2, y = No_of_cells)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group_2, color = Age_group_2)) + 
  stat_summary(aes(x = Age_group_2, y = No_of_cells), fun = median,  fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 30000), breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000)) +
  stat_compare_means(aes(label = paste0("t-test p = ", ..p.adj..), method = "t-test"), 
                     label.y = c(28000), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group") + 
  ylab("Normalize p21 positive cells per 10000 cells per mouse") +  
  ggtitle("p21 positive cells amount per total cells vs age group") +
  scale_color_manual(values=c("#65EFFF", "#996035", "#331900")) + 
  scale_size_manual(values=c(3,3,3,3,3,3)) + 
  theme(text = element_text(size = 20))
