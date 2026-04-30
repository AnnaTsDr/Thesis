library(ggpubr)
library(ggplot2)

TMS_metadata <- read.csv("metedata_TMS.csv")

m <- TMS_metadata %>% group_by(Age_group_2, mouse.id) %>% summarise(No_of_cells = n())

m$Age_group_2 <- factor(m$Age_group_2, levels = c("Young", "Old"))
m$No_of_cells <- as.numeric(as.vector(m$No_of_cells))
my_comparisons <- list( c("Young", "Old"))

ggplot(m, aes(x = Age_group_2, y = No_of_cells)) + 
  geom_jitter(width = 0.15, aes(size = Age_group_2, color = Age_group_2)) + 
  stat_summary(aes(x = Age_group_2, y = No_of_cells), fun = median,  fun.min = median, fun.max = median, 
               geom = "crossbar",
               width = 0.5) + scale_y_continuous(limits = c(0, 30000), breaks = c(0, 5000, 10000, 15000, 20000, 25000,
                                                                                  30000)) +
  stat_compare_means(aes(label = paste0("t-Test p = ", ..p.adj..), method = "t.test"), comparisons = my_comparisons, 
                     label.y = c(27000), label.x = 1.25, size = 5) + 
  theme_classic()+ scale_shape_manual(values=1:23) + xlab("Age group 2") + 
  ylab("Total cells per mouse") +  
  ggtitle("Total cells amount per total cells vs age group 2") +
  scale_color_manual(values=c("#65EFFF", "#996035")) + 
  scale_size_manual(values=c(3,3,3)) + 
  theme(text = element_text(size = 20))


s <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id, subsets)
s <- s %>% summarise(No_of_cells = n())
a <- CD4_T_metadata %>% group_by(Age_group_2, mouse.id) %>% summarise(No_of_total_CD4_T_cells = n())
#normalized total cells number per mice
c1 <- merge(a,s, all = TRUE) %>% filter(subsets == "Naive")
c1$Naive <- "0"
for (n in 1:length(c1[,5])) {
  c1$Naive[n] <- c1$No_of_cells[n]/c1$No_of_total_CD4_T_cells[n]*100
}
c1$Naive <- as.numeric(as.vector(c1$Naive))
c1$subsets <- NULL
c1$No_of_cells <- NULL
c1$No_of_total_CD4_T_cells <- NULL
c2 <- merge(s,a, all = TRUE) %>% filter(subsets == "Naive_Isg15")
c2$Naive_Isg15 <- "0"
for (n in 1:length(c2[,5])) {
  c2$Naive_Isg15[n] <- c2$No_of_cells[n]/c2$No_of_total_CD4_T_cells[n]*100
}
c2$Naive_Isg15 <- as.numeric(as.vector(c2$Naive_Isg15))
c2$subsets <- NULL
c2$No_of_cells <- NULL
c2$No_of_total_CD4_T_cells <- NULL
c <- merge(c1,c2, by = c("mouse.id", "Age_group_2"), all = TRUE)
c3 <- merge(s,a, all = TRUE) %>% filter(subsets == "rTregs")
c3$rTregs <- "0"
for (n in 1:length(c3[,5])) {
  c3$rTregs[n] <- c3$No_of_cells[n]/c3$No_of_total_CD4_T_cells[n]*100
}
c3$rTregs <- as.numeric(as.vector(c3$rTregs))
c3$subsets <- NULL
c3$No_of_cells <- NULL
c3$No_of_total_CD4_T_cells <- NULL
c <- merge(c,c3, by = c("mouse.id", "Age_group_2"), all = TRUE)
c4 <- merge(s,a, all = TRUE) %>% filter(subsets == "aTregs")
c4$aTregs <- "0"
for (n in 1:length(c4[,5])) {
  c4$aTregs[n] <- c4$No_of_cells[n]/c4$No_of_total_CD4_T_cells[n]*100
}
c4$aTregs <- as.numeric(as.vector(c4$aTregs))
c4$subsets <- NULL
c4$No_of_cells <- NULL
c4$No_of_total_CD4_T_cells <- NULL
c <- merge(c,c4, by = c("mouse.id", "Age_group_2"), all = TRUE)
c5 <- merge(s,a, all = TRUE) %>% filter(subsets == "TEM")
c5$TEM <- "0"
for (n in 1:length(c5[,5])) {
  c5$TEM[n] <- c5$No_of_cells[n]/c5$No_of_total_CD4_T_cells[n]*100
}
c5$TEM <- as.numeric(as.vector(c5$TEM))
c5$subsets <- NULL
c5$No_of_cells <- NULL
c5$No_of_total_CD4_T_cells <- NULL
c <- merge(c,c5, by = c("mouse.id", "Age_group_2"), all = TRUE)
c6 <- merge(s,a) %>% filter(subsets == "Exhausted")
c6$Exhausted <- "0"
for (n in 1:length(c6[,5])) {
  c6$Exhausted[n] <- c6$No_of_cells[n]/c6$No_of_total_CD4_T_cells[n]*100
}
c6$Exhausted <- as.numeric(as.vector(c6$Exhausted))
c6$subsets <- NULL
c6$No_of_cells <- NULL
c6$No_of_total_CD4_T_cells <- NULL
c <- merge(c,c6, by = c("mouse.id", "Age_group_2"), all = TRUE)
c7 <- merge(s,a) %>% filter(subsets == "Cytotoxic")
c7$Cytotoxic <- "0"
for (n in 1:length(c7[,5])) {
  c7$Cytotoxic[n] <- c7$No_of_cells[n]/c7$No_of_total_CD4_T_cells[n]*100
}
c7$Cytotoxic <- as.numeric(as.vector(c7$Cytotoxic))
c7$subsets <- NULL
c7$No_of_cells <- NULL
c7$No_of_total_CD4_T_cells <- NULL
c <- merge(c,c7, by = c("mouse.id", "Age_group_2"), all = TRUE)

sen <- TMS_metadata %>% group_by(Age_group_2, mouse.id)
sen <- sen %>% filter(Senescence.group == "Senescence.pos")
sen <- sen %>% summarise(No_of_cells = n())
c8 <- TMS_metadata %>% group_by(Age_group, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
c8 <- merge(sen,c8, all = TRUE)
c8$Senescence <- "0"
for (n in 1:length(c8[,5])) {
  c8$Senescence[n] <- c8$No_of_cells[n]/c8$total_cells[n]*10000
}
c8$Senescence <- as.numeric(as.vector(c8$Senescence))
c8$Senescence.pos <- NULL
c8$No_of_cells <- NULL
c8$total_cells <- NULL
c <- merge(c,c8, by = c("mouse.id", "Age_group_2"), all = TRUE)
c[is.na(c)] <- 0

write.csv(c, "table_T_subsets_senescent_cells.csv")
c$Age_group <- NULL
c <- c %>% filter(mouse.id != "3-M-5/6" & mouse.id != "3-M-7/8" & mouse.id != "3-M-9" & mouse.id != "18-M-52" & 
                    mouse.id != "18-M-53" & mouse.id != "24-M-61" & 
                    mouse.id != "30-M-2" & mouse.id != "3-M-8/9")
