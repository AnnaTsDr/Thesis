
m <- TMS_metadata %>% group_by(Age_group, mouse.id)
m <- m %>% filter(Senescence.group == "Senescence.pos")
m <- m %>% summarise(No_of_senescent_cells = n())
a <- TMS_metadata %>% group_by(Age_group, mouse.id) %>% summarise(total_cells = n())
#normalized total cells number per mice
a <- merge(m,a, all = TRUE)
a[is.na(a)] <- 0
a$normalize <- "0"
for (n in 1:length(a[,5])) {
  a$normalize[n] <- a$No_of_senescent_cells[n]/a$total_cells[n]*10000
}