### EM algorithm for a mixture model of two univariate Gaussian distribution ###

# Load data --------------------------------------------------------------------
data_dir <- "./datasets/em/"

data <- read.table(file.path(data_dir,"Picture.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Picture1.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Picture2.data"),header=TRUE) 
data <- read.table(file.path(data_dir,"Acces.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Barca.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Batman.data"),header=TRUE)
data <- read.table(file.path(data_dir,"BMW.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Face.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Ferrari.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Homer1.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Info.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Lion.data"),header=TRUE)
data <- read.table(file.path(data_dir,"OK.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Porqui.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Robot.data"),header=TRUE)
data <- read.table(file.path(data_dir,"Youtube.data"),header=TRUE)

# show beggining of the data as read # show the pre-parsed data
head(data)

# rename and show the dataset
X<-data$X
X

# some data information
mean(X) # the mean, shows the central expectation of the whole dataset distribution
var(X) # the variance, measure of how much value is away from the mean value.

n<-length(X) # the size, how many individuals?
n

class(X) # what kind of values? numerical variables

### this dataset is an image, we must transform the data as a matrix where each individual corresponds to a pixel
### transform X variable in a matrix, you can rotate image by changing byrow=T o byrow=F
X_image<-matrix(X, nrow=sqrt(n), ncol=sqrt(n), byrow=F) # transpose the unidimensional data array to a bidimensional matrix
class(X_image) # matrix array

# display the matrix numerical values in an image
# the image is just an interpretation of the pixel values
image(X_image)

### we observe we have pixel values
### and we also realize about the distribution of these pixel values observing how the image looks, cuz there are clear differences/variations
### apparently we can assume inside this dataset several sub-population or groups of similars

### observe histogram of the values, to visualize the distribution of the pixel values
### it has a bell shape, composed by random variables
### we may assume the dataset values follows a normal distribution
hist(X_image)


# now, the first thing to do is to guess the values for our initial estimators



# MLE (Maximum Likelihood Estimators)
# One step in the EM algorithm is to obtain MLE of a given statistical model
# The MLE method is a statistical method for estimating the parameters of a statistical model given observations
# This method maximizes the likelihood of making the observations given the parameters of a statistical model

### Start with some initial values for Mu1, Mu2, V1, V2, pi1, pi2
## We could have also considered starting with new values for the Z1, Z2 new random variables, but I prefer to give values to Mu1, Mu2, etc.
Mu1<-185 # expected value for the first sub-population
Mu2<-165 # expected value for the second sub-population
V1<-272.0 # variance of the initial whole dataset
V2<-272.0
pi1<-0.7 # mixture model proportion for the first sub-population
pi2<-1-pi1
####

## Working variables
values<-1
break_s<-0 
pi0<-10.0
n_simul<-1000 ### numbMu1<-round((t(Z1) %*% X)/sum(Z1),4) # max number of iterations if non epsilon converge case
mx1<-matrix(rep(1,n),nrow=n, ncol=1)
epsi<-0.0001 # estimators difference, says when the algorithm converges
######


#Plot to show the evolution of estimates of the parameter pi1
plot(0,pi1,ylim=c(0,1),xlim=c(0,200)) ## in this plot will observe the pi converging, which represents the proportion of belonging to one or another sub-population


###########################LOOP############################
for(j in 1:n_simul) { #First For loop
  if(break_s==0){
    if(values==0){
      ##MLE estimators for a Gaussian mixture model
      Mu1<-round((t(Z1) %*% X)/sum(Z1),4)
      Mu2<-round((t(Z2) %*% X)/sum(Z2),4)
      B1<-X-mx1 %*% Mu1
      B2<-X-mx1 %*% Mu2
      V1<- round((t(B1) %*% diag(array(Z1)) %*% B1)/sum(Z1),4)
      V2<- round((t(B2) %*% diag(array(Z2)) %*% B2)/sum(Z2),4)
      pi1<-round(sum(Z1)/n,4)
      pi2<-round(sum(Z2)/n,4)
    }
    
    ##We always start here for the first iteration
    values<-0 ##parameter values controls that for the first iteration we start here the code
     
    ### calcul new Z1, values of new variables ###
    Exp1<-c() ##Here we set vector Exp1, expectation of the Z1 values
    for(i in 1:n){
      f1<-dnorm(X[i],mean=Mu1,sd=sqrt(V1)) #density function for a Normal variable
      f2<-dnorm(X[i],mean=Mu2,sd=sqrt(V2))
      Exp1[i]<-(pi1*f1)/(pi1*f1+(1.0-pi1)*f2) #Expectation of Z1 values
    }
    ###we can show the Exp1 values in an histogram
    #hist(Exp1)
  
    ### calcul new Z2, values of new variables ###
    ##Note that Exp2 can be just computed as Exp2=1-Exp1
    Exp2<-c()
    for(i in 1:n){
      f1<-dnorm(X[i],mean=Mu1,sd=sqrt(V1))
      f2<-dnorm(X[i],mean=Mu2,sd=sqrt(V2))
      Exp2[i]<-(pi2*f2)/(pi2*f2+(1.0-pi2)*f1) #Expectation of Z2 values
    }

    ####now rewrite Z1 and Z2, these are the values that the new used for new iteration, to obtain new mu1, mu2, V1, v2 and pi1, pi2###
    Z1<-round(Exp1, 4)
    Z2<-round(Exp2, 4)
    print(pi1)

    ##loop to check if we can stop the iteration procedure
    if(abs(pi1-pi0)<epsi){
      break_s<-1
    }
  } else if(break_s==1) {
    break 
  }
  
  pi0<-pi1 
  print(pi1)
  points(j,pi1,pch=19)
}# close First For loop

### Print estimated parameters ###
Mu1; Mu2; V1; V2; pi1

####### Genereta new image ####
### We can display Z1 or Z2 ####

# Z values are the probabilities of belong to one or another sub-population

X_image2<-matrix(t(Z1),nrow=sqrt(n),ncol=sqrt(n),byrow=F)
image(X_image2)

### assign 0 or 1 in terms of Z1 probabilities ###
Z11<-array(rep(0,n))
for(i in 1:n){
  if(Z1[i]>=pi1){
    Z11[i]=1
  }
}
Z11_image<-matrix(Z11,nrow=sqrt(n),ncol=sqrt(n),byrow=F)

image(Z11_image)


### To generate a graphic ###
postscript(file="Image3.ps")
layout(matrix(c(1,2),1,2,byrow=TRUE),respect=FALSE)
image(X_image)
image(Z11_image)
dev.off()
