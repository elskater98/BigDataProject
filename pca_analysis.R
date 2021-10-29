library(dplyr)
library(ggplot2)
library(corrplot)
library(factoextra)

# Load data
data_dir <- "./datasets/pca/"

# Read data
data = read.csv(file.path(data_dir, 'BodyFat.csv'), sep = ',')

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
R <- cor(centrate_data, method="pearson")
corrplot(R, method="color")

### Compute the S matrix (variance-covariance matrix of X) ###
# variance-covariance matrix of X
S <- (t(centrate_data) %*% centrate_data) / (nrow(data))

##Remenber the variance-covariance matrix are the poblational parameters (/n), try sum(B[,1]^2)/n=S[1,1], sum(B[,2]^2)/n=S[2,2], or
##sum(B[,1]*B[,2])/n=S[1,2]=S[2,1]

### Obtain a diagonal matrix D_1s of order m with the inver of standard deviation (using S)
D_1s <- diag(sqrt(1 / diag(S)), nrow = ncol(data), ncol = ncol(data))

### Now D_1s is a diagonal matrix of order m with the inver of standard deviation
### Obtain the Y matrix (centrered and reduced X matrix) ###
Y <- centrate_data %*% D_1s # Centered matrix and reduced

###Note that variance-covariance matrix of Y is equal to the R Matrix
###Try sum(Y[,1]^2)/n=R[1,1], sum(Y[,2]^2)/n=R[2,2], or sum(Y[,1]*Y[,2])/n=R[1,2]=R[2,1]

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

#ggplot(data = data.frame(prop_varianza_acum, pc = 1:ncol(data)),
#       aes(x = pc, y = prop_varianza_acum, group = 1)) +
#  geom_point() +
#  geom_line() +
#  theme_bw() +
#  labs(x = "Componente principal",
#       y = "Prop. varianza explicada acumulada")

ggplot(data = data.frame(prop_varianza_acum, pc = factor(1:ncol(data))),
       aes(x = pc, y = prop_varianza_acum, group = 1)) +
  geom_point() +
  geom_line() +
  geom_label(aes(label = round(prop_varianza_acum,2))) +
  theme_bw() +
  labs(x = "Componentes principales", 
       y = "Prop. varianza explicada acumulada")

fviz_eig(pca)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)