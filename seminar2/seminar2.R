rnorm(n = 10)
rnorm(10)

set.seed(540)
rnorm(10)

n <- 10
B <- 4

set.seed(540)
x <- matrix(rnorm(n * B), nrow = n)
str(x)
head(x)

set.seed(540)
x <- matrix(rnorm(n * B), nrow = n)
rownames(x) <- sprintf("obs%02d", 1:n)
colnames(x) <- sprintf("samp%02d", 1:B)
x
dim(x)
mean(x[,2])
colMeans(x)
apply(x, 2, mean)

# compute the average of the sample means
mean(colMeans(x))

# try with more observations and samples
n = 100
B = 10
set.seed(540)
x <- matrix(rnorm(n * B), nrow = n)
colMeans(x)
mean(colMeans(x))