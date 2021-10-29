### EM algorithm for a mixture model of two uni-variate Gaussian distribution ###

# Load data
data_dir <- "./datasets/em/"

## Read data
data <- read.table(file.path(data_dir,"Picture.data"))
data <- read.table(file.path(data_dir,"Picture1.data"))
data <- read.table(file.path(data_dir,"Picture2.data"))
data <- read.table(file.path(data_dir,"Acces.data"))
data <- read.table(file.path(data_dir,"Barca.data"))
data <- read.table(file.path(data_dir,"Batman.data"))
data <- read.table(file.path(data_dir,"BMW.data"))
data <- read.table(file.path(data_dir,"Face.data"))
data <- read.table(file.path(data_dir,"Ferrari.data"))
data <- read.table(file.path(data_dir,"Homer1.data"))
data <- read.table(file.path(data_dir,"Info.data"))
data <- read.table(file.path(data_dir,"Lion.data"))
data <- read.table(file.path(data_dir,"OK.data"))
data <- read.table(file.path(data_dir,"Porqui.data"))
data <- read.table(file.path(data_dir,"Robot.data"))
data <- read.table(file.path(data_dir,"Youtube.data"))

# show beginning of the data as read # show the pre-parsed data
head(data)

# show the structure of the data (type, shape, number of variables, kind of variables)
str(data)

X<-data$X # getting and assigning the variable X from the data.frame to work directly with an array of data
X

### the mean and variance are very important
### they provide a global view of the values
mean(X) # the mean, shows the central expectation of the whole data-set distribution, among all possible inner clusters
var(X) # the variance, measure of how much value is away from the mean/expectation value. the range.

n<-length(X) # the size, how many individuals?
n

summary(X)

class(X) # what kind of values? numerical variables

### an image is basically a large interpreted amount of numbers
### this data-sets aims to represent an image, where each individual corresponds to a pixel value
### commonly there are RGB images based on three channel pixel representation
### but here, we have only one channel as there is only one variable, therefore, one value for pixel

### this data-set is an image, we must transform the data as a matrix where each individual corresponds to a pixel
### transform X variable in a matrix, you can rotate image by changing byrow=T o byrow=F
X_image<-matrix(X, nrow=sqrt(n), ncol=sqrt(n), byrow=F) # transpose the uni-dimensional data array to a bi-dimensional matrix
class(X_image) # type: matrix array
X_image

# display the matrix numerical values in an image
# the image is just an interpretation of the pixel values
image(X_image) # represent the matrix array into an image 

### we observe we have pixel values
### the color of the pixel is proportional to the value of that pixel
### image interpretation => darker color means higher value, lighter color means lower value

### also, we realize about the distribution of these pixel values observing how the image looks, because there are clear differences/variations
### apparently we can assume inside this data-set several sub-population or groups of similars

### it's possible to misclassify some pixel values, as there are pixel exceptions into the color areas recognized
### sometimes is so difficult to differentiate exceptions and assign them to the correct cluster
### but this how works statistics, as we are working with probabilities
### for instance, we can't say all women and all men have distinct height => this is not true, as maybe there are similars, becoming undifferentiated in values terms

### objective: classify the image pixel values
### what kind of analysis may we apply? The EM Analysis

### one important thing to display when trying to classify, is the histogram of the data-set
### observe the histogram of the values, to visualize the distribution/range of the pixel values in numerical/statistical terms
### it has a bell shape, composed by random variables
### we may assume the data-set values follows a normal distribution
hist(X_image)

### from the histogram we may assume there are sub-populations inside the data-set
### observing the distribution, the symmetry, and the variance of the values
### in some case will be apparent the clusters, but others, may be not so clear

### now, the first thing to do for the EM Analysis is to guess/give the values for our initial estimators based on our pre-analysis/observations
### that's the reason why the histogram is useful, as observing the distribution in numerical values terms

# MLE (Maximum Likelihood Estimators)
# One step in the EM algorithm is to obtain MLE of a given statistical model
# The MLE method is a statistical method for estimating the parameters of a statistical model given observations
# This method maximizes the likelihood of making the observations given the parameters of a statistical model

### depending on how many cluster we identified, we will assume a mixture model
### from that one, we will need more or less estimators variables

### as we identified two clusters,
### start with some initial values for Mu1, Mu2, V1, V2, pi1, pi2 for the first iteration assumption
### we could have also considered starting with new values for the Z1, Z2 new random variables, but I prefer to give values to Mu1, Mu2, etc. Easier
Mu1<-185 # expected value for the first sub-population
Mu2<-165 # expected value for the second sub-population
V1<-272# variance of the initial whole data-set
V2<-272
pi1<-0.7 # mixture model proportion for the first sub-population
pi2<-1-pi1

### now let's try if the EM algorithm may converge
### thus, the EM algorithm will have classified the pixel values
### and we will able to see if the classification obtained is similar two the pre-analysis observed
### also, must be interesting how this initial values may affect the resulting classification

# working variables
values<-1 # first iteration control
break_s<-0 # epsilon converge control
pi0<-10.0
n_simul<-1000 ### numbMu1<-round((t(Z1) %*% X)/sum(Z1),4) # max number of iterations if non epsilon converge case
mx1<-matrix(rep(1,n),nrow=n, ncol=1)
epsi<-0.0001 # estimators difference, says when the algorithm converges, establishes the permissive differentiation

# Plot to show the evolution of estimates of the parameter pi1
plot(x=0, y=pi1, ylim=c(0,1), xlim=c(0,200)) ## in this plot will observe the pi converging, which represents the proportion of belonging to one or another sub-population

###########################LOOP############################
for(j in 1:n_simul) { # First For loop: controls de maximum number of EM iterations allowed
  if(break_s==0) {
    if(values==0) { # from the second iteration on, we must approximate the estimators values based on the previous Expectation
      
      ### now rewrite Z1 and Z2, these are the values that the new used for new iteration, to obtain new mu1, mu2, V1, v2 and pi1, pi2
      ### the new Z values are based on the previous iteration Expectation calculated
      ### Z(1) => Mu1(1) => Exp(1) => Z(2) => Mu1(2) => Exp(2) => ... etc.
      Z1<-round(Exp1, 4) # Z values are probabilistic indirect proportional variables that classify our data-set into two clusters
      Z2<-round(Exp2, 4)
      
      ## MLE estimators for a Gaussian mixture model
      Mu1<-round((t(Z1) %*% X)/sum(Z1),4) # maximization step
      Mu2<-round((t(Z2) %*% X)/sum(Z2),4) # maximization step
      B1<-X-mx1 %*% Mu1
      B2<-X-mx1 %*% Mu2
      V1<-round((t(B1) %*% diag(array(Z1)) %*% B1)/sum(Z1),4)
      V2<-round((t(B2) %*% diag(array(Z2)) %*% B2)/sum(Z2),4)
      pi1<-round(sum(Z1)/n, 4) # calculate the probability proportion
      pi2<-round(sum(Z2)/n, 4)
    }
    else {
      values<-0 ## controls that for the first iteration we start here the code
    }
    
    ## We always start here for the first iteration
    ## we calculate for the first iteration the Z values and the expectation for each sub-population, based on our pre-defined estimators values
    ## the ones are coherent with our pre-analysis when observing the data-set distribution
     
    ### calculate new Z1, values of new variables ###
    Exp1<-c() ## Here we set vector Exp1, expectation of the Z1 values
    for(i in 1:n) {
      f1<-dnorm(X[i],mean=Mu1,sd=sqrt(V1)) # density function for a Normal variable
      f2<-dnorm(X[i],mean=Mu2,sd=sqrt(V2))
      Exp1[i]<-(pi1*f1)/(pi1*f1+(1.0-pi1)*f2) # Expectation of Z1 values
    }
    
    ### we can show the Exp1 values in an histogram
    #hist(Exp1)
  
    ### calculate new Z2, values of new variables
    Exp2<-c()
    for(i in 1:n){
      f1<-dnorm(X[i],mean=Mu1,sd=sqrt(V1))
      f2<-dnorm(X[i],mean=Mu2,sd=sqrt(V2))
      Exp2[i]<-(pi2*f2)/(pi2*f2+(1.0-pi2)*f1) # Expectation of Z2 values
    }
    
    ## note that Exp2 can be just computed as Exp2=1-Exp1
    #Exp2=1-Exp1
    
    # plot update the evolution of pi (proportion of the mixture model)
    points(j, pi1, pch=19) # x = iteration number, y = absolute pi value, pch = specify the plotting character
    print(paste("pi1:", pi1))
    
    ## loop to check if we can stop the iteration procedure
    if(abs(pi1-pi0) < epsi){
      break_s<-1
      print("The algorithm has converged!!")
    }
    
    # copy the new pi value obtained to be used for checking the next iteration distance differentiation (epsilon)
    pi0<-pi1
    
  } else if(break_s==1) {
    break 
  }
}# close First For loop

### Print estimated parameters ###
Mu1; Mu2; V1; V2; pi1;

hist(Z1)
hist(Z2)

####### Generate new image ####
### We can display Z1 or Z2 ####
# Z values are the probabilities of belong to one or another sub-population

X_image2<-matrix(t(Z1),nrow=sqrt(n),ncol=sqrt(n),byrow=F)
image(X_image2)

### assign 0 or 1 in terms of Z1 probabilities ###
### checking results
Z11<-array(rep(0,n))
data2<-read.table(file.path(data_dir,"Batman.data"))
X2<-data2$X

for(i in 1:n){
  if(Z1[i]>=pi1) {
    Z11[i]=1
  }
}

Z11_image<-matrix(Z11,nrow=sqrt(n),ncol=sqrt(n),byrow=F)
image(Z11_image)

### To generate a graphic ###
postscript(file="Image3.ps")
layout(matrix(c(1,2),1,2, byrow=TRUE), respect=FALSE)
image(X_image)
image(Z11_image)
dev.off()