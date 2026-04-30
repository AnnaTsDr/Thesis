library(ggpubr)
library(ggplot2)

#Average subsets pie
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
summarise(a, mean = mean())

mouse_50 <- a %>% filter(mouse.id == "18-F-50")
pie(as.vector(mouse_50$percent), labels = mouse_50$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_50&midpoint <- cumsum(percent) - percent / 2
p <- ggplot(mouse_50, aes(x="", y=percent, fill=subsets))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + theme_minimal()+ theme(
    axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(),
    panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"),
    axis.text.x=element_blank())
p + geom_text_repel(aes(label = percent))


mouse_51 <- a %>% filter(mouse.id == "18-F-51")
pie(as.vector(mouse_51$percent), labels = mouse_51$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_52 <- a %>% filter(mouse.id == "18-M-52")
pie(as.vector(mouse_52$percent), labels = mouse_52$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_53 <- a %>% filter(mouse.id == "18-M-53")
pie(as.vector(mouse_53$percent), labels = mouse_53$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_54 <- a %>% filter(mouse.id == "21-F-54")
pie(as.vector(mouse_54$percent), labels = mouse_54$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_55 <- a %>% filter(mouse.id == "21-F-55")
pie(as.vector(mouse_55$percent), labels = mouse_55$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_58 <- a %>% filter(mouse.id == "24-M-58")
pie(as.vector(mouse_58$percent), labels = mouse_58$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_59 <- a %>% filter(mouse.id == "24-M-59")
pie(as.vector(mouse_59$percent), labels = mouse_59$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_60 <- a %>% filter(mouse.id == "24-M-60")
pie(as.vector(mouse_60$percent), labels = mouse_60$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_61 <- a %>% filter(mouse.id == "24-M-61")
pie(as.vector(mouse_61$percent), labels = mouse_61$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_2 <- a %>% filter(mouse.id == "30-M-2")
pie(as.vector(mouse_2$percent), labels = mouse_2$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_3 <- a %>% filter(mouse.id == "30-M-3")
pie(as.vector(mouse_3$percent), labels = mouse_3$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_4 <- a %>% filter(mouse.id == "30-M-4")
pie(as.vector(mouse_4$percent), labels = mouse_4$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_5 <- a %>% filter(mouse.id == "30-M-5")
pie(as.vector(mouse_5$percent), labels = mouse_5$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_62 <- a %>% filter(mouse.id == "1-M-62")
pie(as.vector(mouse_62$percent), labels = mouse_62$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_63 <- a %>% filter(mouse.id == "1-M-62")
pie(as.vector(mouse_63$percent), labels = mouse_63$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_56 <- a %>% filter(mouse.id == "3-F-56")
pie(as.vector(mouse_56$percent), labels = mouse_56$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_57 <- a %>% filter(mouse.id == "3-F-57")
pie(as.vector(mouse_57$percent), labels = mouse_57$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_5_6 <- a %>% filter(mouse.id == "3-M-5/6")
pie(as.vector(mouse_5_6$percent), labels = mouse_5_6$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_7_8 <- a %>% filter(mouse.id == "3-M-7/8")
pie(as.vector(mouse_7_8$percent), labels = mouse_7_8$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_8 <- a %>% filter(mouse.id == "3-M-8")
pie(as.vector(mouse_8$percent), labels = mouse_8$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))
mouse_9 <- a %>% filter(mouse.id == "3-M-9")
pie(as.vector(mouse_9$percent), labels = mouse_9$subsets, 
    col = c("red","green","blue","yellow","grey","pink","orange"))





ggplot(a, aes(x = Age_group, y = normalize)) + 
  geom_jitter(width = 0.15, aes(shape = mouse.id, size = Age_group)) + stat_summary(aes(x = Age_group, y = normalize), fun = median,  fun.min = median, fun.max = median, geom = "crossbar", 
                                                                                    width = 0.5) + scale_y_continuous(limits = c(0, 630), breaks = c(0, 630)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(570, 600, 630)) + 
  theme_classic()+ scale_shape_manual(values=1:22)

write.csv(a, "CD4_for_pie_average.csv")
