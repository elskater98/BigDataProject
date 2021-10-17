library(dplyr)
library(ggplot2)
library(corrplot)

data = read.csv('datasets/pca/School_qualifications2.csv', sep = ',')

# Data Information
names(data) # Column names
sapply(data, typeof) # Column types
dim(data) # Dimensions
summary(data) # Summary
View(data) # View data as table



# Step by step PCA
# https://davetang.org/muse/2012/02/01/step-by-step-principal-components-analysis-using-r/
# https://rpubs.com/Joaquin_AR/287787
# Correlation Covariance, Variance: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
centrate_data <- data

# Mean for each column (Standardization)
for (i in names(centrate_data)) {
  centrate_data[i] <- data[[i]] - mean(data[[i]])
  
}

centrate_data <- as.matrix(centrate_data) # dataframe to matrix

# Correlation Matrix
R <- cor(centrate_data, method = "pearson")
corrplot(R, method="color")

# variance-covariance matrix of X
S <- (t(centrate_data) %*% centrate_data) / (nrow(data))

D_1s <-
  diag(sqrt(1 / diag(S)), nrow = ncol(data), ncol = ncol(data))

Y <- centrate_data %*% D_1s # Centered matrix and reduced

### Obtain Eigenvalues and Eigenvector from matrix R ###

Eig <- eigen(R)

### Eigenvalues ###

Eig$values

### Eigenvectors ###

Eig$vectors

# R Library
pca <- prcomp(centrate_data, center = TRUE, scale. = TRUE)
pca$rotation # Eigen Vectors
pca$sdev ^ 2 #Eigen Values
pca$x # Principal Components
prop_varianza <- pca$sdev^2/sum(pca$sdev^2)
prop_varianza_acum <- cumsum(prop_varianza)

biplot(x = pca, scale = 0, cex = 0.6, col = c("blue4", "brown3"))

ggplot(data = data.frame(prop_varianza, pc = 1:4),
       aes(x = pc, y = prop_varianza)) +
  geom_col(width = 0.3) +
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  labs(x = "Componente principal",
       y = "Prop. de varianza explicada")

ggplot(data = data.frame(prop_varianza_acum, pc = 1:ncol(data)),
       aes(x = pc, y = prop_varianza_acum, group = 1)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "Componente principal",
       y = "Prop. varianza explicada acumulada")

ggplot(data = data.frame(prop_varianza_acum, pc = factor(1:ncol(data))),
       aes(x = pc, y = prop_varianza_acum, group = 1)) +
  geom_point() +
  geom_line() +
  geom_label(aes(label = round(prop_varianza_acum,2))) +
  theme_bw() +
  labs(x = "Componentes principales", 
       y = "Prop. varianza explicada acumulada")
