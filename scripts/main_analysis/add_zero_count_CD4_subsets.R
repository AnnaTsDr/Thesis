CD4_T_metadata <- read.csv("CD4_T_metedata_TMS.csv")
m <- CD4_T_metadata %>% group_by(Age_group, mouse.id, subsets)
m <- m %>% summarise(No_of_cells = n())
m <- as.data.frame(m)
m <- m %>% complete(mouse.id, subsets)
m <- m %>% group_by(Age_group, mouse.id, subsets)
for (i in 1:154) {
  ifelse(grep("^[13]-", m[i,1]), m[i,3] <- "Young", m[i,3] <- "Old")
  ifelse(grep("[8240]-", m[i,1]), m[i,3] <- "Old", m[i,3] <- "Young")
}
m[is.na(m)] <- 0
m <- m[,c(3, 1, 2, 4)]
write.csv(m, "CD4_T_subsets.csv")

TMS_metadata <- read.csv("metedata_TMS.csv")
m <- TMS_metadata %>% group_by(Age_group, mouse.id, Senescence.group)
m <- m %>% summarise(No_of_cells = n())
m <- as.data.frame(m)
m <- m %>% tidyr::complete(mouse.id, Senescence.group)
m <- m %>% group_by(Age_group, mouse.id, Senescence.group)
for (i in 1:45) {
  ifelse(grep("Senescence.pos", m[i,1]), m[i,3] <- "Young", m[i,3] <- "Old")
  ifelse(grep("Senescence.*", m[i,1]), m[i,3] <- "Old", m[i,3] <- "Young")
}
m[is.na(m)] <- 0
m <- m[,c(3, 1, 2, 4)]
write.csv(m, "CD4_T_subsets.csv")