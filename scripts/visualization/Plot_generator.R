mdata <- read.csv("metadata_TMS_clistered.csv")

ggplot(mdata, aes(sex, Senescence.group)) + geom_boxplot()

m <- mdata %>% group_by(Age_group, mouse.id)
m <- m %>% summarise()
m <- m %>% summarise(No_of_mice = n())
a <- mdata %>% group_by(Age_group, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a)
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$total_cells[n]/a$No_of_mice[n]*10
}
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old", "supercentenerians"))
a$normalize <- as.numeric(as.vector(a$normalize))
my_comparisons <- list( c("Young", "Old"), c("Young", "supercentenerians"), c("Old", "supercentenerians"))
ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 64000), breaks = c(0, 64000)) +
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(50000, 55000, 60000)) + 
  theme_classic()
#there is no significant difference between the age groups

#not normalized total cells
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old", "supercentenerians"))
a$total_cells <- as.numeric(as.vector(a$total_cells))
my_comparisons <- list( c("Young", "Old"), c("Young", "supercentenerians"), c("Old", "supercentenerians"))
ggplot(a, aes(x = Age_group, y = total_cells, size = 1, color = Age_group)) + 
  geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group, y = total_cells), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 26000), breaks = c(0, 26000)) +
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(20000, 22000, 24000)) + 
  theme_classic()
#there is no significant difference between the age groups

#young density plot
y <- a %>% filter(Age_group == "Young")
ggplot(y, aes(x = total_cells)) + geom_density() + theme_classic() + 

#old density plot
o <- a %>% filter(Age_group == "Old")
ggplot(o, aes(x = total_cells)) + geom_density() + theme_classic()

#supercentenerians density plot
s <- a %>% filter(Age_group == "supercentenerians")
ggplot(s, aes(x = total_cells)) + geom_density() + theme_classic()

#3m 18m 24m 30m
age_3_18_24_30 <- mdata %>% group_by(Age_group, age, mouse.id) %>% summarise(total_cells = n()) %>% filter(age == "3m" | age == "18m" | age == "24m" | 
                                                                                                             age == "30m")
ggplot(age_3_18_24_30, aes(x = total_cells, color = age)) + geom_density() + theme_classic()

descdist(y, discrete = FALSE)#to idetify the destribution

#T cells of TMS
a$Age_group[11:14] <- "supercentenarians"
m_t <- read.csv("metadata_t_cells.csv")
t_cells <- m_t %>% group_by(Age_group, mouse.id, my_clusters) %>% summarise(No_type_cells = n()) %>% filter(my_clusters == "CD8 T" | 
                                                                                                                my_clusters == "CD4 T" |
                                                                                                              my_clusters == "DP T" |
                                                                                                              my_clusters == "DN T")
cd4 <- t_cells %>% filter(my_clusters == "CD4 T")
cd8 <- t_cells %>% filter(my_clusters == "CD8 T")
cd4_cd8 <- merge(cd4, cd8, by = "mouse.id")
t_cells <- merge(a, t_cells)
t_cells$cd4 <- cd4
t_cells <- merge(t_cells, cd8, all = TRUE)
cd4_cd8$cd4_cd8 <- "0"
for (n in 1:length(cd4_cd8[,5])) {
  cd4_cd8$cd4_cd8[n] <- cd4_cd8$No_type_cells.x[n]/cd4_cd8$No_type_cells.y[n]*100
}
cd4_cd8$Age_group.x <- factor(cd4_cd8$Age_group.x, levels = c("Young", "Old", "supercentenarians"))
cd4_cd8$cd4_cd8 <- as.numeric(as.vector(cd4_cd8$cd4_cd8))
my_comparisons <- list( c("Young", "Old"), c("Young", "supercentenarians"), c("Old", "supercentenarians"))
ggplot(cd4_cd8, aes(x = Age_group.x, y = cd4_cd8)) + 
  geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group.x, y = cd4_cd8), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 400), breaks = c(0, 400)) +
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(320, 360, 400)) + 
  theme_classic()
ggplot(cd4_cd8, aes(x = Age_group.x, y = No_type_cells.x)) + 
  geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group.x, y = No_type_cells.x), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 800), breaks = c(0, 800)) +
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(720, 760, 800)) + 
  theme_classic()

cd4 <- mdata %>% group_by(Age_group, mouse.id, free.annotation) %>% filter(free.annotation == "CD4+ T") %>% summarise(No_type_cells = n())
cd8 <- mdata %>% group_by(Age_group, mouse.id, free.annotation) %>% filter(free.annotation == "CD8+ T") %>% summarise(No_type_cells = n())
cd4 %>% rename(CD4 = No_type_cells)
cd8 %>% rename(CD8 = No_type_cells)
cd4_cd8 <- merge(a, cd4)
cd4_cd8 <- merge(cd4, cd8)
cd4_cd8$cd4_cd8 <- "0"
for (n in 1:length(cd4_cd8[,5])) {
  cd4_cd8$cd4_cd8[n] <- cd4_cd8$CD4[n]/cd4_cd8$CD8[n]*100
}
ggplot(cd4_cd8, aes(x = Age_group, y = cd4_cd8)) + 
  geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group, y = cd4_cd8), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 800), breaks = c(0, 800)) +
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(720, 760, 800)) + 
  theme_classic()



t_cells$normalized <- "0"
for (n in 1:length(t_cells[,5])) {
  t_cells$normalized[n] <- t_cells$No_type_cells[n]/t_cells$total_cells[n]*10000
}
t_cells$Age_group <- factor(t_cells$Age_group, levels = c("Young", "Old", "supercentenarians"))
t_cells$normalized <- as.numeric(as.vector(t_cells$normalized))
my_comparisons <- list( c("Young", "Old"), c("Young", "supercentenarians"), c("Old", "supercentenarians"))
ggplot(t_cells, aes(x = Age_group, y = normalized)) + 
  geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group, y = normalized), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                           width = 0.5) + scale_y_continuous(limits = c(0, 1010), breaks = c(0, 1010)) +
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(900, 950, 1000)) + 
  theme_classic()
t_cells$cd4_cd8 <- "0"
for (n in 1:length(t_cells[,5])) {
  t_cells$cd4_cd8[n] <- t_cells$cd4[n]/t_cells$cd8[n]*100
}





m_m <- read.csv("t_or_macrophages_subset_metadata.csv")
mph_cells <- m_m %>% 



b <- mdata %>% group_by(tissue, Age_group, mouse.id, Senescence.group) %>% summarise(n = n())

y <- a %>% filter(Age_group == "Young")

ggplot(a, aes(n)) + geom_histogram(binwidth = 100) + stat_bin(bins = 1)
ggplot(y, aes(n)) + geom_histogram(binwidth = 100) + stat_bin(bins = 3)
mdata %>% group_by(Age_group, mouse.id) %>% filter(Age_group == "Young") %>% summarise(n = n()) %>% group_by(Age_group) %>% 
  summarise(max = max(n), min = min(n), mean = mean(n), median = median(n)) %>% ggplot(aes(Age_group)) + geom_boxplot() 
mdata %>% group_by(Age_group, mouse.id, Senescence.group, T_cells, Monocytes) %>% filter(Age_group == "Young") %>% stat_boxplot(aes(Age_group))

FeatureScatter()
ggplot(b, aes(x = mouse.id, group = Age_group)) + geom_jitter(aes(y = median)) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
                                                                                              geom = "crossbar", width = 0.5) + 
  geom_dotplot(aes(x = Senescence.group))

a <- b %>% filter(Senescence.group == "Senescence.pos")
a$Age_group <- factor(a$Age_group, levels = c("Young", "Old", "supercentenerians"))#reorder x axis
ggplot(a, aes(x = Age_group, y = n)) + geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group, y = n),fun = median, fun.min = median, fun.max = median, 
                                                                                              geom = "crossbar", width = 0.5) 


# Statistical test
stat.test <- compare_means(n ~ Age_group, method = "wilcox.test", data = a)
stat.test

c <- mdata %>% group_by(tissue, Age_group, mouse.id) %>% summarise(total_cells = n())
d <- merge(a, c, all = TRUE)
d$normalised <- "0"
d[is.na(d)] <- "0"
for (n in 1:length(d[,6])) {
  d$normalized[n] <- d$n[n]/d$total_cells[n]*10000
}
d$Age_group <- factor(d$Age_group, levels = c("Young", "Old", "supercentenerians"))#reorder x axis
ggplot(d, aes(x = Age_group, y = total_cells)) + geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group, y = total_cells), fun = median, fun.min = median, fun.max = median, 
                                                                                geom = "crossbar", width = 0.5) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(20000, 23000, 26000)) + theme_classic()


# Statistical test
stat.test <- compare_means(total_cells ~ Age_group, method = "t.test", data = d)
stat.test

d$nornalized <- as.numeric(as.vector(d$normalized))
my_comparisons <- list( c("Young", "Old"), c("Young", "supercentenerians"), c("Old", "supercentenarians"))
ggplot(d, aes(x = Age_group, y = normalized, group = Age_group)) + geom_jitter(width = 0.15) + stat_summary(aes(x = Age_group, y = normalized, group = Age_group), fun = median, fun.min = median, fun.max = median, 
                                                                                          geom = "crossbar", width = 0.5) + 
scale_y_continuous(limits = c(0, 200), breaks = c(0, 50, 100, 150, 200)) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(100, 130, 160)) + theme_classic()
    
# Statistical test
stat.test <- compare_means(normalized ~ Age_group, method = "t.test", data = d)
stat.test

#Normalized senescent cells per mice per age (1m, 3m, 18m, 21m, 24m, 30m)
md_group_by_age <- mdata %>% group_by(age, mouse.id, Senescence.group) %>% summarise(senecent_cells = n()) %>% filter(Senescence.group == "Senescence.pos")
md_group_by_age_total <- mdata %>% group_by(age, mouse.id) %>% summarise(total_cells = n())
md_group_by_age <- merge(md_group_by_age, md_group_by_age_total, all = TRUE)
md_group_by_age$normalized <- "0"
md_group_by_age[is.na(md_group_by_age)] <- "0"
for (n in 1:length(md_group_by_age[,6])) {
  md_group_by_age$normalized[n] <- md_group_by_age$senecent_cells[n]/md_group_by_age$total_cells[n]*10000
}
md_group_by_age$normalized <- as.numeric(as.vector(md_group_by_age$normalized))
md_group_by_age$age <- factor(md_group_by_age$age, levels = c("1m", "3m", "18m","21m", "24m", "30m"))#reorder x axis
ggplot(md_group_by_age, aes(x = age, y = normalized, group = age)) + geom_jitter(width = 0.15) + stat_summary(aes(x = age, y = normalized, group = age), fun = median, fun.min = median, fun.max = median, 
                                                                                                     geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 200), breaks = c(0, 50, 100, 150, 200)) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = list(c("1m", "3m"), c("1m", "18m"), c("1m", "21m"), c("1m", "24m"),
                                                                                      c("1m", "30m"), c("3m", "18m"), c("3m", "21m"), c("3m", "24m"), 
                                                                                      c("3m", "30m"), c("18m", "21m"), c("18m", "24m"), c("18m", "30m"),
                                                                                      c("21m", "24m"), c("21m", "30m"), c("24m", "30m")), 
                     label.y = c(60, 70, 80 ,90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)) + theme_classic()

# Statistical test
stat.test <- compare_means(normalized ~ age, method = "t.test", data = md_group_by_age)
stat.test

#
c <- mdata %>% group_by(tissue, Age_group, mouse.id) %>% summarise(total_cells = n())
a <- mdata %>% group_by(tissue, Age_group, mouse.id, Senescence.group) %>% summarise(senescent_cells = n()) %>% filter(Senescence.group == "Senescence.pos")
md_group_by_tissue <- merge(a, c, all = TRUE) %>% filter(Senescence.group == "Senescence.pos")
md_group_by_tissue$normalized <- "0"
md_group_by_tissue[is.na(md_group_by_tissue)] <- "0"
for (n in 1:length(md_group_by_tissue[,6])) {
  md_group_by_tissue$normalized[n] <- md_group_by_tissue$senescent_cells[n]/md_group_by_tissue$total_cells[n]*10000
}
md_group_by_tissue$normalized <- as.numeric(as.vector(md_group_by_tissue$normalized))
md_group_by_tissue$Age_group <- factor(md_group_by_tissue$Age_group, levels = c("Young", "Old", "supercentenerians"))
ggplot(md_group_by_tissue, aes(x = Age_group, y = normalized, group = Age_group)) + 
  geom_jitter(width = 0.15) + 
  stat_summary(aes(x = Age_group, y = normalized, group = Age_group), fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 850), breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850)) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, 
                     label.y = c(500, 600, 700)) + theme_classic()
stat.test <- compare_means(normalized ~ Age_group, method = "t.test", data = md_group_by_tissue)



c <- mdata %>% group_by(tissue, Age_group, mouse.id) %>% summarise(total_cells = n())
a <- mdata %>% group_by(tissue, Age_group, mouse.id, Senescence.group) %>% summarise(senescent_cells = n()) %>% filter(Senescence.group == "Senescence.pos")
md_group_by_tissue <- merge(a, c, all = TRUE) %>% filter(Senescence.group == "Senescence.pos")
md_group_by_tissue$normalized <- "0"
md_group_by_tissue[is.na(md_group_by_tissue)] <- "0"
for (n in 1:length(md_group_by_tissue[,6])) {
  md_group_by_tissue$normalized[n] <- md_group_by_tissue$senescent_cells[n]/md_group_by_tissue$total_cells[n]*10000
}
md_group_by_tissue$normalized <- as.numeric(as.vector(md_group_by_tissue$normalized))
md_group_by_tissue$tissue <- factor(md_group_by_tissue$tissue, levels = c("Bladder", "Fat", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", 
                                                                          "Lung", "Spleen", "Marrow", "Pancreas", "Thymus", "Large_Intestine", 
                                                                          "Mammary_Gland", "Skin", "Tongue", "Trachea"))
md_group_by_tissue <- md_group_by_tissue %>% filter(Age_group == "Young")
ggplot(md_group_by_tissue, aes(x = tissue, y = normalized, group = tissue)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = tissue, y = normalized, group = tissue), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 300), breaks = c(0, 50, 100, 150, 200, 250, 300)) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", 
                     label.y = c(100, 130, 160)) + theme_classic()
stat.test <- compare_means(normalized ~ tissue, method = "t.test", data = md_group_by_tissue)


md_group_by_tissue <- merge(a, c, all = TRUE)  %>% filter(Senescence.group == "Senescence.pos")
md_group_by_tissue$normalized <- "0"
md_group_by_tissue[is.na(md_group_by_tissue)] <- "0"
for (n in 1:length(md_group_by_tissue[,6])) {
  md_group_by_tissue$normalized[n] <- md_group_by_tissue$senescent_cells[n]/md_group_by_tissue$total_cells[n]*10000
}
md_group_by_tissue$normalized <- as.numeric(as.vector(md_group_by_tissue$normalized))
md_group_by_tissue$tissue <- factor(md_group_by_tissue$tissue, levels = c("Bladder", "Fat", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", 
                                                                          "Lung", "Spleen", "Marrow", "Pancreas", "Thymus", "Large_Intestine", 
                                                                          "Mammary_Gland", "Skin", "Tongue", "Trachea"))
md_group_by_tissue <- md_group_by_tissue %>% filter(Age_group == "Old")
ggplot(md_group_by_tissue, aes(x = tissue, y = normalized, group = tissue)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = tissue, y = normalized, group = tissue), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 850), breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850)) + theme_classic()# + 
#  stat_compare_means(aes(label = ..p.signif..), method = "t.test", 
   #                  label.y = c(100, 130, 160))
#stat.test <- compare_means(normalized ~ tissue, method = "t.test", data = md_group_by_tissue)


md_group_by_tissue <- merge(a, c, all = TRUE)
md_group_by_tissue$normalized <- "0"
md_group_by_tissue[is.na(md_group_by_tissue)] <- "0"
for (n in 1:length(md_group_by_tissue[,6])) {
  md_group_by_tissue$normalized[n] <- md_group_by_tissue$senescent_cells[n]/md_group_by_tissue$total_cells[n]*10000
}
md_group_by_tissue$normalized <- as.numeric(as.vector(md_group_by_tissue$normalized))
md_group_by_tissue$tissue <- factor(md_group_by_tissue$tissue, levels = c("Bladder", "Fat", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", 
                                                                          "Lung", "Spleen", "Marrow", "Pancreas", "Thymus", "Large_Intestine", 
                                                                          "Mammary_Gland", "Skin", "Tongue", "Trachea"))
md_group_by_tissue <- md_group_by_tissue %>% filter(Age_group == "supercentenerians")
ggplot(md_group_by_tissue, aes(x = tissue, y = normalized, group = tissue)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = tissue, y = normalized, group = tissue), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 300), breaks = c(0, 50, 100, 150, 200, 250, 300)) + theme_classic()
 # stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, 
  #                   label.y = c(100, 130, 160)) + 
stat.test <- compare_means(normalized ~ tissue, method = "t.test", data = md_group_by_tissue)

md_group_by_tissue <- merge(a, c, all = TRUE) %>% filter(tissue == "Kidney")
md_group_by_tissue$normalized <- "0"
for (n in 1:length(md_group_by_tissue[,6])) {
  md_group_by_tissue$normalized[n] <- md_group_by_tissue$senescent_cells[n]/md_group_by_tissue$total_cells[n]*10000
}
md_group_by_tissue[is.na(md_group_by_tissue)] <- "0"
md_group_by_tissue$normalized <- as.numeric(as.vector(md_group_by_tissue$normalized))
md_group_by_tissue$tissue <- factor(md_group_by_tissue$tissue, levels = c("Bladder", "Fat", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", 
                                                                          "Lung", "Spleen", "Marrow", "Pancreas", "Thymus", "Large_Intestine", 
                                                                          "Mammary_Gland", "Skin", "Tongue", "Trachea"))
md_group_by_tissue$Age_group <- factor(md_group_by_tissue$Age_group, levels = c("Young", "Old", "supercentenerians"))
ggplot(md_group_by_tissue, aes(x = Age_group, y = normalized, group = Age_group)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = Age_group, y = normalized, group = Age_group), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) + theme_classic() + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(60, 80, 100)) 
stat.test <- compare_means(normalized ~ tissue, method = "t.test", data = md_group_by_tissue)



c <- mdata %>% group_by(tissue, mouse.id) %>% summarise(total_cells = n())
c$tissue <- factor(c$tissue, levels = c("Bladder", "Fat", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", 
                                                                          "Lung", "Spleen", "Marrow", "Pancreas", "Thymus", "Large_Intestine", 
                                                                          "Mammary_Gland", "Skin", "Tongue", "Trachea"))
ggplot(c, aes(x = tissue, y = total_cells, group = tissue)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = tissue, y = total_cells, group = tissue), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 10000), breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", 
                     label.y = c(100, 130, 160)) + theme_classic()
stat.test <- compare_means(total_cells ~ tissue, method = "t.test", data = c)



#Marrow analysis
md_group_by_tissue <- merge(a, c, all = TRUE) %>% filter(tissue == "Marrow")
md_group_by_tissue$normalized <- "0"
for (n in 1:length(md_group_by_tissue[,6])) {
  md_group_by_tissue$normalized[n] <- md_group_by_tissue$senescent_cells[n]/md_group_by_tissue$total_cells[n]*10000
}
#md_group_by_tissue[is.na(md_group_by_tissue)] <- "0"
md_group_by_tissue$normalized <- as.numeric(as.vector(md_group_by_tissue$normalized))
md_group_by_tissue$Age_group <- factor(md_group_by_tissue$Age_group, levels = c("Young", "Old", "supercentenerians"))
ggplot(md_group_by_tissue, aes(x = Age_group, y = normalized, group = Age_group)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = Age_group, y = normalized, group = Age_group), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) + theme_classic() + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(60, 80, 100)) 
stat.test <- compare_means(normalized ~ tissue, method = "t.test", data = md_group_by_tissue)

Marrow <- readRDS("F:/Immunaging/Clustred_Marrow.rds")
Marrow[["Age_group"]] <- plyr::mapvalues(x = Marrow$age, from = c("1m", "3m", "18m", "21m", "24m", "30m"),
                                            to = c("Young", "Young", "Old", "Old", "Old", "supercentenarians"))
Marrow_md <- Marrow@meta.data
Marrow_md <- Marrow_md %>% group_by(Age_group, mouse.id) %>% summarise(total_cells = n())
Marrow_T <- Marrow@meta.data %>% group_by(Age_group, mouse.id, type) %>% summarise(T_cells = n()) %>% filter(type == "T cell-Marrow")
Marrow_md_T <- merge(Marrow_md, Marrow_T, all = TRUE) %>% filter(type == "T cell-Marrow")
Marrow_md_T$normalized <- "0"
for (n in 1:length(Marrow_md_T[,6])) {
  Marrow_md_T$normalized[n] <- Marrow_md_T$T_cells[n]/Marrow_md_T$total_cells[n]*10000
}
Marrow_md_T$normalized <- as.numeric(as.vector(Marrow_md_T$normalized))
Marrow_md_T$Age_group <- factor(Marrow_md_T$Age_group, levels = c("Young", "Old", "supercentenarians"))
ggplot(Marrow_md_T, aes(x = Age_group, y = normalized, group = Age_group)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = Age_group, y = normalized, group = Age_group), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 600), breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600)) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, 
                     label.y = c(100, 130, 160)) + theme_classic()
stat.test <- compare_means(normalized ~ Age_group, method = "t.test", data = Marrow_md_T)


Marrow_M <- Marrow@meta.data %>% group_by(Age_group, mouse.id, type) %>% summarise(Monocyte_cells = n()) %>% filter(type == "monocyte-Marrow" | 
                                                                                                               type == "macrophage-Marrow")
Marrow_md_M <- merge(Marrow_md, Marrow_M, all = TRUE) %>% filter(type == "monocyte-Marrow" | type == "macrophage-Marrow")
Marrow_md_M$normalized <- "0"
for (n in 1:length(Marrow_md_M[,6])) {
  Marrow_md_M$normalized[n] <- Marrow_md_M$Monocyte_cells[n]/Marrow_md_M$total_cells[n]*10000
}
Marrow_md_M$normalized <- as.numeric(as.vector(Marrow_md_M$normalized))
Marrow_md_M$Age_group <- factor(Marrow_md_M$Age_group, levels = c("Young", "Old", "supercentenarians"))
ggplot(Marrow_md_M, aes(x = Age_group, y = normalized, group = Age_group)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = Age_group, y = normalized, group = Age_group), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 1300), breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300)) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, 
                     label.y = c(900, 1000, 1100)) + theme_classic()
stat.test <- compare_means(normalized ~ Age_group, method = "t.test", data = Marrow_md_M)

#Spleen
md_group_by_tissue <- merge(a, c, all = TRUE) %>% filter(tissue == "Spleen")
md_group_by_tissue$normalized <- "0"
for (n in 1:length(md_group_by_tissue[,6])) {
  md_group_by_tissue$normalized[n] <- md_group_by_tissue$senescent_cells[n]/md_group_by_tissue$total_cells[n]*10000
}
#md_group_by_tissue[is.na(md_group_by_tissue)] <- "0"
md_group_by_tissue$normalized <- as.numeric(as.vector(md_group_by_tissue$normalized))
md_group_by_tissue$Age_group <- factor(md_group_by_tissue$Age_group, levels = c("Young", "Old", "supercentenerians"))
ggplot(md_group_by_tissue, aes(x = Age_group, y = normalized, group = Age_group)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = Age_group, y = normalized, group = Age_group), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) + theme_classic() + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(60, 80, 100)) 
stat.test <- compare_means(normalized ~ tissue, method = "t.test", data = md_group_by_tissue)

Spleen <- readRDS("F:/Immunaging/Clustred_Spleen.rds")
Spleen[["Age_group"]] <- plyr::mapvalues(x = Spleen$age, from = c("1m", "3m", "18m", "21m", "24m", "30m"),
                                         to = c("Young", "Young", "Old", "Old", "Old", "supercentenarians"))
Spleen_md <- Spleen@meta.data
Spleen_md <- Spleen_md %>% group_by(Age_group, mouse.id) %>% summarise(total_cells = n())
Spleen_T <- Spleen@meta.data %>% group_by(Age_group, mouse.id, type) %>% summarise(T_cells = n()) %>% filter(type == "T cell-Spleen")
c <- c %>% filter(tissue == "Spleen")
Spleen_md_T <- merge(Spleen_md, Spleen_T) 
Spleen_md_T$normalized <- "0"
for (n in 1:length(Spleen_md_T[,6])) {
  Spleen_md_T$normalized[n] <- Spleen_md_T$T_cells[n]/Spleen_md_T$total_cells[n]*10000
}
Spleen_md_T$normalized <- as.numeric(as.vector(Spleen_md_T$normalized))
Spleen_md_T$Age_group.y <- factor(Spleen_md_T$Age_group, levels = c("Young", "Old", "supercentenarians"))
ggplot(Spleen_md_T, aes(x = Age_group.y, y = normalized, group = Age_group.y)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = Age_group.y, y = normalized, group = Age_group.y), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 3500), breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500)) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, 
                     label.y = c(60000, 65000, 70000)) + theme_classic()
stat.test <- compare_means(normalized ~ Age_group.y, method = "t.test", data = Spleen_md_T)


Spleen_M <- Spleen@meta.data %>% group_by(Age_group, mouse.id, type) %>% summarise(Monocyte_cells = n()) %>% filter(type == "monocyte-Marrow" | 
                                                                                                                      type == "macrophage-Spleen")
c <- c %>% filter(type == "monocyte-Marrow" | type == "macrophage-Spleen") 
Spleen_md_M <- merge(Spleen_md, Spleen_M)# %>% filter(type == "monocyte-Marrow" | type == "macrophage-Spleen")
Spleen_md_M$normalized <- "0"
for (n in 1:length(Spleen_md_M[,6])) {
  Spleen_md_M$normalized[n] <- Spleen_md_M$Monocyte_cells[n]/Spleen_md_M$total_cells[n]*10000
}
Spleen_md_M$normalized <- as.numeric(as.vector(Spleen_md_M$normalized))
Spleen_md_M$Age_group.y <- factor(Spleen_md_M$Age_group, levels = c("Young", "Old", "supercentenarians"))
ggplot(Spleen_md_M, aes(x = Age_group.y, y = normalized, group = Age_group.y)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = Age_group.y, y = normalized, group = Age_group.y), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 1500), breaks = c(0, 1500)) + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, 
                     label.y = c(900, 1000, 1100)) + theme_classic()
stat.test <- compare_means(normalized ~ Age_group, method = "t.test", data = Spleen_md_M)

#Spleen after reclustering





#Liver
md_group_by_tissue <- merge(a, c, all = TRUE) %>% filter(tissue == "Liver")
md_group_by_tissue$normalized <- "0"
for (n in 1:length(md_group_by_tissue[,6])) {
  md_group_by_tissue$normalized[n] <- md_group_by_tissue$senescent_cells[n]/md_group_by_tissue$total_cells[n]*10000
}
md_group_by_tissue[is.na(md_group_by_tissue)] <- "0"
md_group_by_tissue$normalized <- as.numeric(as.vector(md_group_by_tissue$normalized))
md_group_by_tissue$tissue <- factor(md_group_by_tissue$tissue, levels = c("Bladder", "Fat", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", 
                                                                          "Lung", "Spleen", "Marrow", "Pancreas", "Thymus", "Large_Intestine", 
                                                                          "Mammary_Gland", "Skin", "Tongue", "Trachea"))
md_group_by_tissue$Age_group <- factor(md_group_by_tissue$Age_group, levels = c("Young", "Old", "supercentenerians"))
ggplot(md_group_by_tissue, aes(x = Age_group, y = normalized, group = Age_group)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = Age_group, y = normalized, group = Age_group), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 250), breaks = c(0, 50, 100,150, 200, 250)) + theme_classic() + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(60, 80, 100)) 
stat.test <- compare_means(normalized ~ tissue, method = "t.test", data = md_group_by_tissue)

#Lung
md_group_by_tissue <- merge(a, c, all = TRUE) %>% filter(tissue == "Lung")
md_group_by_tissue$normalized <- "0"
for (n in 1:length(md_group_by_tissue[,6])) {
  md_group_by_tissue$normalized[n] <- md_group_by_tissue$senescent_cells[n]/md_group_by_tissue$total_cells[n]*10000
}
md_group_by_tissue[is.na(md_group_by_tissue)] <- "0"
md_group_by_tissue$normalized <- as.numeric(as.vector(md_group_by_tissue$normalized))
md_group_by_tissue$tissue <- factor(md_group_by_tissue$tissue, levels = c("Bladder", "Fat", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", 
                                                                          "Lung", "Spleen", "Marrow", "Pancreas", "Thymus", "Large_Intestine", 
                                                                          "Mammary_Gland", "Skin", "Tongue", "Trachea"))
md_group_by_tissue$Age_group <- factor(md_group_by_tissue$Age_group, levels = c("Young", "Old", "supercentenerians"))
ggplot(md_group_by_tissue, aes(x = Age_group, y = normalized, group = Age_group)) + 
  geom_jitter(width = 0.15, mapping = aes(color = mouse.id)) + 
  stat_summary(aes(x = Age_group, y = normalized, group = Age_group), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5) + 
  scale_y_continuous(limits = c(0, 150), breaks = c(0, 50, 100, 150)) + theme_classic() + 
  stat_compare_means(aes(label = ..p.signif..), method = "t.test", comparisons = my_comparisons, label.y = c(60, 80, 100)) 
stat.test <- compare_means(normalized ~ tissue, method = "t.test", data = md_group_by_tissue)
