data = read.csv('datasets/pca/School_qualifications2.csv',sep = ',')
#as.matrix(data)

# Data Information
names(data) # Column names
sapply(data, typeof) # Column types
dim(data) # Dimensions
summary(data) # Summary
View(data) # View data as table


# Step by step PCA
# https://davetang.org/muse/2012/02/01/step-by-step-principal-components-analysis-using-r/
# https://rpubs.com/Joaquin_AR/287787

centrate_data <- data

# Mean for each column (Standardization)
for (i in names(centrate_data)) {
  centrate_data[i] <- data[[i]] - mean(data[[i]])
} 


View(centrate_data) # B

# Covariance Matrix (Correlation)
matrix_cov <- cov(centrate_data)

View(matrix_cov)

# Eigen Vectors and Values
eigen <- eigen(matrix_cov)

# Transpose Data
t_eigenvectors <- t(eigen$vectors)
t_centrate_data <- t(centrate_data)

# Matrix Product
pc_scores <- t_eigenvectors %*% t_centrate_data
rownames(pc_scores) <- names(data)
View(t(pc_scores))

# PCA using R Library
pca <- prcomp(centrate_data, scale = TRUE, center = TRUE)
View(pca$scale)
pca$sdev

biplot(x = pca, scale = 0, cex = 0.6, col = c("blue4", "brown3"))

prop_varianza <- pca$sdev^2/sum(pca$sdev^2)
prop_varianza_acum <- cumsum(prop_varianza)


