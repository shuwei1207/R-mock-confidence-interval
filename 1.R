devtools::install_git("https://github.com/ccolonescu/PoEdata")
library(PoEdata)
data("andy")
mod2 <- lm(sales~price+advert+I(advert^2),data=andy)
smod2 <- summary(mod2)

b1 <- coef(mod2)[[1]]
b2 <- coef(mod2)[[2]]
b3 <- coef(mod2)[[3]]
b4 <- coef(mod2)[[4]]


//第一部分觀察資料
alpha <- 0.05
df <- mod2$df.residual
tcr <- qt(1-alpha/2, df)
g <- (1-b3)/(2*b4)
g3 <- -1/(2*b4)
g4 <- -(1-b3)/(2*b4^2)
varb3 <- vcov(mod2)[3,3]
varb4 <- vcov(mod2)[4,4]
covb3b4 <- vcov(mod2)[3,4]
varg <- g3^2*varb3+g4^2*varb4+2*g3*g4*covb3b4
seg <- sqrt(varg)
lowbg <- g-tcr*seg
upbg <- g+tcr*seg


//第二部分從模型生成資料

N<-40
alpha <- 0.05
df <- mod2$df.residual
nrsim<-1000
x2<-andy$price
x3<-andy$advert
tcr <- qt(1-alpha/2, df)
cri<-(1-b3)/(2*b4)

vb1<-numeric(nrsim)
vb2<-numeric(nrsim)
vb3<-numeric(nrsim)
vb4<-numeric(nrsim)

g<-numeric(nrsim)
g3<-numeric(nrsim)
g4<-numeric(nrsim)
varb3<-numeric(nrsim)
varb4<-numeric(nrsim)
covb3b4<-numeric(nrsim)
varg<-numeric(nrsim)
seg<-numeric(nrsim)
lowbg<-numeric(nrsim)
upbg<-numeric(nrsim)


for(i in 1:nrsim){
	set.seed(12345+10*i)
	y <- b1+b2*x2+b3*x3+b4*x3*x3+rnorm(N,mean=0,sd=1)
	mod3 <- lm(y~x2+x3+I(x3^2))
	vb1[i] <- coef(mod3)[[1]]
	vb2[i] <- coef(mod3)[[2]]
	vb3[i] <- coef(mod3)[[3]]
	vb4[i] <- coef(mod3)[[4]]
	g[i]<-(1-vb3[i])/(2*vb4[i])
	g3[i] <- -1/(2*vb4[i])
	g4[i] <- -(1-vb3[i])/(2*vb4[i]^2)
	varb3[i] <- vcov(mod3)[3,3]
	varb4[i] <- vcov(mod3)[4,4]
	covb3b4[i] <- vcov(mod3)[3,4]
	varg[i] <- g3[i]^2*varb3[i]+g4[i]^2*varb4[i]+2*g3[i]*g4[i]*covb3b4[i]
	seg[i] <- sqrt(varg[i])
	lowbg[i] <- g[i]-tcr*seg[i]
	upbg[i] <- g[i]+tcr*seg[i]
}

//第三部分資料若落在區間則=1 資料若不在則=0 求總和
a<-numeric(nrsim)
for(i in 1:nrsim)
{
	if(lowbg[i]<cri&&cri<upbg[i])a[i]=1
	else a[i]=0
}
result=sum(a)


//第四部份結論
