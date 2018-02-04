# Test matrix subtraction
nr <- 10000
nc <- 30
m <- matrix(sample(nr*nc), nrow=nr)

repeats <- 1000
subtracted.row <- m[1, ]
t0 <- proc.time()
for (i in 1:repeats) { rowSums(sweep(m, 2, subtracted.row)^2) }
t1 <- proc.time()
for (i in 1:repeats) { rowSums((m - rep(subtracted.row, each=nr))^2) }
t2 <- proc.time()

print(t1-t0)
print(t2-t1)
