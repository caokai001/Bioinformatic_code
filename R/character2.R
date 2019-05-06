###2
x<-c(1,3,2,5)
x
x=c(1,6,2)
x
y=c(1,4,3)
length(x)
length(y)
x+y
ls()
rm(x,y)
ls()
rm(list=ls())
?matrix
x=matrix(data=c(1,2,3,4),nrow=2,ncol=2)
x
x=matrix(c(1,2,3,4),2,2)
x
sqrt(x)
x^2
x=rnorm(50)
y=x+rnorm(50,mean=50,sd=0.1)
cor(x,y)
set.seed(1303)
rnorm(50)
set.seed(3)
y=rnorm(100)
mean(y)
var(y)
sqrt(var(y))
sd(y)
###2.3.2
x=rnorm(100)
y=rnorm(100)
plot(x,y)
plot(x,y,xlab="x ",ylab="y",main="plot of x and y")
pdf("Figure.pdf")
plot(x,y,col="green")
dev.off()
x=seq(1,10)
x
x=1:10
x=seq(-pi,pi,length=50)
x
y=x
f=outer(x,y,function(x,y)cos(y)/(1+x^2))
contour(x,y,f)
contour(x,y,f,nlevels = 45,add=T)
fa=(f-t(f))/2
contour(x,y,fa,nlevels = 15)
image(x,y,fa)
persp(x,y,fa)
persp(x,y,fa,theta = 30,phi=40)
###2.3.3
A=matrix(1:16,4,4)
A
A[2,3]
A[,c(1,2)]
A[1,]
A[-c(1,3)]
dim(A)
###2.3.4
library(ISLR)



#####3. linear Regression
###3.6.1
library(MASS)
library(ISLR)
fix(Bosten)
names(Boston)
lm.fit=lm(medv~lstat,data=Boston)
attach(Boston)
lm.fit=lm(medv~lstat)
lm.fit
summary(lm.fit)
names(lm.fit)
coef(lm.fit)
confint(lm.fit)
predict(lm.fit,data.frame(lstat=(c(5,10,15))),interval ="confidence" )
plot(lstat,medv)
abline(lm.fit,col="red")
plot(1:20,1:20,pch=1:20)
par(mfrow=c(2,2))
plot(lm.fit)
plot(predict(lm.fit),residuals(lm.fit))
plot(predict(lm.fit),rstudent(lm.fit))
###3.6.3
lm.fit=lm(medv~lstat+age,data=Boston)
summary(lm.fit)
lm.fit=lm(medv~.,data=Boston)
summary(lm.fit)
library(car)
vif(lm.fit)   #方差膨胀因子
lm.fit1=lm(medv~.-age,data=Boston)
summary(lm.fit1)
lm.fit1=update(lm.fit,~.-age)
###3.6.4
summary(lm(medv~lstat*age,data=Boston))
###3.6.5
lm.fit2=lm(medv~lstat+I(lstat^2))
summary(lm.fit2)
lm.fit=lm(medv~lstat)

anova(lm.fit,lm.fit2)
par(mfrow=c(2,2))
plot(lm.fit2)
lm.fit5=lm(medv~poly(lstat,5))
summary(lm.fit5)
summary(lm(medv~log(rm),data=Boston))
####3.6.6
fix(Carseats)
names(Carseats)
lm.fit=lm(Sales~.+Income:Advertising+Price:Age,data=Carseats)
summary(lm.fit)
attach(Carseats)
contrasts(ShelveLoc)
###3.6.7
LoadLibraries=function(){
  library(ISLR)
  library(MASS)
  print("The libraries have been loaded.")
}
LoadLibraries()
