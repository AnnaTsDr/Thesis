p_a <- list()
p_b <- list()

x <- 1
for (i in 1:17) {
  for (j in i+1:17) {
    a <- rownames(cd4t0n6.markers[[i]])
    b <- rownames(cd4t0n6.markers[[j]])
    c <- intersect(a,b)
    p_a[[x]] <- (length(c)/length(a))*100
    p_b[[x]] <- (length(c)/length(b))*100
    if (p_a[[x]] > 50 | p_b[[x]] > 50){ mark [i,1] <- i
    mark [i,2] <- j}
    x <- x+1
  }
}

a <- rownames(cd4t0n6.markers_1)
b <- rownames(cd4t0n6.markers_3)
c <- intersect(a,b)
p_a <- (length(c)/length(a))*100
p_b <- (length(c)/length(b))*100