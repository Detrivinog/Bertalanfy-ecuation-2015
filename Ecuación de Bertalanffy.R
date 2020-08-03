setwd("C:/Users/David Esteban/Google Drive/Materias cursadas/Ecología Forestal/Taller")
datos<-read.csv2("datos.csv")
colnames(datos)=c("D","dD")

with(datos, plot(D,dD))
#Función estimar.m
estimar.m<-function(dw,w,datos,x,y,z){
        m=seq(x,y,z)
        CME<-function(mod){
                cme= CME=anova(mod)$Mean[length(anova(mod)$Mean)]
        }
        cme=vector()
        for(i in 1:length(m)){
                mod<-lm(dw~I(w^m[i])+w-1,datos)
                cme[i]=CME(mod)
        } 
        plot(m,cme,type="l", xlab="m",ylab="CME del modelo lineal")
        data=data.frame(cme,m)
        print(data[which.min(data$cme),]) 
}
#Grafica del resultado de estimar m
estimar.m(datos$dD,datos$D,datos,1.16,1.18,0.0001)
    
estimar.m(datos$dD,datos$D,datos,1.14,1.19,0.001)

with(m, plot(m,cme,xlab="m",ylab="CME del modelo lineal"))

# 0.007389716 1.1708

#Modelo de Bertalanffy dD/dt
m=1.1708
mod<-lm(dD~I(D^m)+D-1,datos)
summary(mod)
anova(mod)
n=summary(mod)$coefficients[1,1] #Catabolismo
g=summary(mod)$coefficients[2,1] #Anabolismo

#Grafica Modelo Bertalanffy dD/dt
op <- par(mar = c(5,4.5,4,2) + 0.1)
with(datos, plot(D,dD,xlab="DAP (cm)", ylab=expression(dDAP/dt~(cm/año)), 
                 xlim=c(0,72), ylim=c(0,1.4) ))
par(new=TRUE)
curve(n*(x^m)+g*x, xlab="",col="red", 
      ylab="", xlim=c(0,72), ylim=c(0,1.4), axes=FALSE)

#Capacidad de carga dD/dt=0
dD=n*(datos$D^m)-g*datos$D
k=(g/(-1*n))^(1/(m-1))

#Tasa de crecimiento en la que se máximiza la función
dD2=(-g/(n*m))^(1/(m-1)) 
dDmax=n*(dD2^m)+g*dD2

#Integrada 
A=((-n/g)^(1/m-1))
b=1-(1.1/A)^(1/(m-1))
t=seq(0,100,1)
D=k*(1+b*exp(g*(1-m)*t))^(1/(1-m))
plot(D, type="l")

D<-function(t, m, n, g ){
        k=(g/(-1*n))^(1/(m-1))
        A=((-n/g)^(1/m-1))
        b=1-(1.1/A)^(1/(m-1))
       D=k*(1+b*exp(g*(1-m)*t))^(1/(1-m))
       print(D)
}

d<-D(t, m, n, g)
plot(d, type="l")
plot(t)
tcarga=log(((((.99*k)/k)^(1-m))-1)/b)/(g*(1-m))

#Modelo relativo
datos$dD_D<-with(datos, dD/D)
mr=1.87639
nr=-0.00140648
gr=0.0641619
modr<-lm(dD_D~I(D^(mr-1)), datos)
summary(modr)

dDr<-function(x){
        mr=1.87639
        -0.00140648*(x^(mr-1))+0.0641619}
with(datos, plot(D,dD_D, xlab="D (cm)", ylab="1/D(dD/dt)",pch="*",
                 xlim=c(0,72), ylim=c(0,0.12)))
par(new=TRUE)
curve(dDr,1,72, xlab="",col="blue", 
      ylab="", xlim=c(0,72), ylim=c(0,0.12), axes=FALSE)

#Residuales
plot(predict(modr),rstudent(modr), ylab="Residuales estudentizados",
     xlab="Predichos", ylim=c(-5,5))
abline(h=0, col="blue")
curve(dDr,1,72)

