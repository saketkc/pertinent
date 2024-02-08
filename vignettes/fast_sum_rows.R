
counts.matrix <- matrix(c(1,0,0,0, 1,0,0,0, 0,1,2,3, 3,1,0,0, 0,0,0,1 ), ncol=4)
counts.matrix <- as(counts.matrix, "dgCMatrix")

# > counts.matrix
#      [,1] [,2] [,3] [,4]
# [1,]    1    0    2    0
# [2,]    0    0    3    0
# [3,]    0    0    3    0
# [4,]    0    0    1    0
# [5,]    1    1    0    1
row.names <- c("chr1", "chr1", "chr2", "chr2", "chr3")
mymodel.matrix <- sparse.model.matrix(
  object = ~ 0 + row.names
)
colnames(x = mymodel.matrix) <- sapply(
  X = colnames(x = mymodel.matrix),
  FUN = function(name) {
    name <- gsub(pattern = "row.names", replacement = "", x = name)
    return(paste0(rev(x = unlist(x = strsplit(x = name, split = ":"))),
                  collapse = "__"
    ))
  }
)
# > mymodel.matrix
# 5 x 3 sparse Matrix of class "dgCMatrix"
# row.nameschr1 row.nameschr2 row.nameschr3
# 1             1             .             .
# 2             1             .             .
# 3             .             1             .
# 4             .             1             .
# 5             .             .             1
rownames(counts.matrix) <- row.names
counts.agg <- as.matrix(x = (t(mymodel.matrix) %*% counts.matrix))
# > counts.agg
# [,1] [,2] [,3] [,4]
# chr1    1    0    5    0
# chr2    0    0    4    0
# chr3    1    1    0    1