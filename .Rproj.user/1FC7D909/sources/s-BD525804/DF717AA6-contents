
############################################
######  FUNCTIONS - life table functions
############################################

###### female
LT_reconstruction_female <- function(mx) {
  age<-(0:(length(mx)-1)) ;  ax<-rep(0.5,length(mx))
  qx<- fx <- lx <- Lx <- Tx <- ex <-c()
  mx[is.na(mx)]=0
  if (mx[1] < 0.01724) {
    ax[1]=0.14903-2.05527*mx[1]
  } else {
    if (mx[1] < 0.06891){
      ax[1]=0.04667 + 3.88089*mx[1]
    } else {
      ax[1]=0.31411
    }
  }
  lastNonNA <- max(which((!is.na(mx)) & (mx!=0) ))
  if (lastNonNA>100){
    ax[lastNonNA]=1/mx[lastNonNA]
  }
  if (lastNonNA<length(age)){
    ax<-ax[-c(lastNonNA+1:length(age))]
    mx<-mx[-c(lastNonNA+1:length(age))]
    age<-age[-c(lastNonNA+1:length(age))]
  }
  qx=(mx)/(1+(1-ax)*mx) ; lx<-c(1,cumprod(1-qx)) ; lx<-lx[1:length(mx)] ; lx[lx<0]=0
  fx=lx*qx ; Lx=lx-fx+ax*fx ; Tx=rev(cumsum(rev(Lx))) ; ex=Tx/lx ; ex[is.na(ex)]=0
  LT<-data.frame(Age=age, ax=ax, mx=mx, qx=qx, fx=fx, lx=lx, Lx=Lx, Tx=Tx, ex=ex)
  return(LT)
}
###### male
LT_reconstruction_male <- function(mx) {
  age<-(0:(length(mx)-1)) ; ax<-rep(0.5,length(mx))
  qx<- fx <- lx <- Lx <- Tx <- ex <-c()
  mx[is.na(mx)]=0
  if (mx[1] < 0.023) {
    ax[1]=0.14929-1.99545*mx[1]
  } else {
    if (mx[1] < 0.08307){
      ax[1]=0.02832 + 3.26021*mx[1]
    } else {
      ax[1]=0.29915
    }
  }
  lastNonNA <- max(which((!is.na(mx)) & (mx!=0) ))
  if (lastNonNA>100){
    ax[lastNonNA]=1/mx[lastNonNA]
  }
  if (lastNonNA<length(age)){
    ax<-ax[-c(lastNonNA+1:length(age))]
    mx<-mx[-c(lastNonNA+1:length(age))]
    age<-age[-c(lastNonNA+1:length(age))]
  }
  qx=(mx)/(1+(1-ax)*mx) ; lx<-c(1,cumprod(1-qx)) ; lx<-lx[1:length(mx)] ; lx[lx<0]=0
  fx=lx*qx ; Lx=lx-fx+ax*fx ; Tx=rev(cumsum(rev(Lx))) ; ex=Tx/lx ; ex[is.na(ex)]=0
  LT<-data.frame(Age=age, ax=ax, mx=mx, qx=qx, fx=fx, lx=lx, Lx=Lx, Tx=Tx, ex=ex)
  return(LT)
}
