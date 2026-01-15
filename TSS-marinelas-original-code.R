library("MASS")
library("glmnet")
library("changepoint")
library("SIS")
library("ncvreg")
library("stabs")
library("car")
library("plsRglm")
library("rJava")
library("glmulti")
library("spikeslab")

    

prostate=as.matrix(read.csv('H:/Biostatistics/Profiles/CapanuM/multiple regression/Large p/HD paper/Data analyses/Prostate/SIS prostate/prostate_2000.csv', sep = ',', header=TRUE))


# define response variable
x=std(prostate[,-1]);y=prostate[,1]


#x is the matrix of predictors
#y is the outcome

cutoff=0.625 #the cutoff used in Stage 1, obtained using the stabsel_parameters function;
#stabsel_parameters(p=2000,q=10,PFER=0.2, assumptions="none")

nrepeat=50 #number of subsamples in Stage 1
nrepeat2=50 #number of subsamples in Stage 2
percent_split=0.5 #take subsamples of half size
nsim=1 #number of times you run the method
family="binomial" #type of outcome



tss=function(x,y,nrepeat,nrepeat2,percent_split,nsim,family,cutoff){

full_data=cbind(y,x)

m=ncol(x)
n=nrow(x)
m_restr=round(n/(4*log(n)))

k_cv=matrix(NA,nrow=nsim,ncol=m)
k_cv_fin2=matrix(NA,nrow=nsim,ncol=m)
signif_threshold=rep(NA,nsim)
signif1=rep(NA,nsim)

TSS_NA=rep(0,nsim);

failTSS=0;

#sink(nullfile()) 

for (isim in 1:nsim) {print(isim);
 
freq_fin=matrix(0,nrow=nrepeat, ncol=m)

freq_union=matrix(0,nrow=nrepeat, ncol=m)



for (jsim in 1:nrepeat) {#print("jsim");print(jsim)

 
  split=sample(1:n,round(percent_split*n),replace=FALSE)
   train_data=full_data[split,]
   test_data=full_data[-split,]
   #first column is y and then followed by the m predictors

   y_train=train_data[,1]
   x_train=train_data[,2:(m+1)]

  
scad_subs=cv.ncvreg(x_train,y_train,family=family,penalty="SCAD")
  fit_scad_subs=scad_subs$fit
  beta_scad_subs=fit_scad_subs$beta[,scad_subs$min]


freq_fin[jsim,which(beta_scad_subs[2:(m+1)]!=0)]=1

   nm=1:m
   left_1=nm[-which(beta_scad_subs[2:(m+1)]!=0)];

    try.pls <- try(
    modplsglm <- plsRglm(y_train,x_train[,left_1],2,modele="pls-glm-logistic")
  )

  if(class(try.pls) != "try-error"){
    
if (!is.null(modplsglm$pp)) {

if (ncol(modplsglm$pp)==2) {
    pc1_order=order(-abs(modplsglm$pp[,1]))
    pc2_order=order(-abs(modplsglm$pp[,2]))
   
 
   index1_sel=pc1_order[1:m_restr];
   index2_sel=pc2_order[1:m_restr];
  
   freq_union[jsim,union(left_1[union(index1_sel,index2_sel)],which(beta_scad_subs[2:(m+1)]!=0))]=1
                           }

    else {if (ncol(modplsglm$pp)==1) {TSS_NA[isim]=1;
          pc1_order=order(-abs(modplsglm$pp[,1]))
          index1_sel=pc1_order[1:m_restr];
  
  
   freq_union[jsim,union(left_1[index1_sel],which(beta_scad_subs[2:(m+1)]!=0))]=1
                           }
         }

                } else {freq_union[jsim,which(beta_scad_subs[2:(m+1)]!=0)]=1; TSS_NA[isim]=1;}

} else {freq_union[jsim,which(beta_scad_subs[2:(m+1)]!=0)]=1; TSS_NA[isim]=1;}

#need a } to close the nrepeat loop is it here?
}



k_cv[isim,]= apply(freq_union,2,mean)


signif=which(k_cv[isim,]>=cutoff)
signif1[isim]=length(signif)

if (length(signif)>0) {

freq_fin2=matrix(0,nrow=nrepeat, ncol=m)
kcv_fin2=rep(NA,nrepeat)



for (jsim2 in 1:nrepeat2) {

 
  split2=sample(1:n,round(percent_split*n),replace=FALSE)
   train_data=full_data[split2,]
   test_data=full_data[-split2,]
   #first column is y and then followed by the m predictors

   y_train=train_data[,1]
   x_train=train_data[,2:(m+1)]

  regrdata2=as.data.frame(cbind(y_train,x_train[,signif]))   

glm2=glm(y_train ~ ., data =regrdata2,family=binomial)
bestgl2=do.call("glmulti", list(glm2, method="h", level=1, crit="aic",plotty=F,report=F,family=family))
bestreg2=bestgl2@objects[[1]]
keepx2=names(bestreg2$model)[-1]
idsel2=match(keepx2, names(regrdata2))-1

freq_fin2[jsim2,signif[idsel2]]=1

#need a } to close the nrepeat loop is it here?
}


k_cv_fin2[isim,]= apply(freq_fin2,2,mean)


sort_decrease=order(k_cv_fin2[isim,],decreasing=TRUE)

try.out <- try(
    changepointTSS <- cpt.meanvar(k_cv_fin2[isim,sort_decrease],method="PELT") 
  )

  if(length(signif)>2 & class(try.out) != "try-error"){
cptsTSS=cpts(changepointTSS);
cpts_TSS=cptsTSS[1]
signif_TSS=which(k_cv_fin2[isim,]>=k_cv_fin2[isim,sort_decrease[cpts_TSS]])
signif_threshold[isim]=k_cv_fin2[isim,sort_decrease[cpts_TSS]]



if (length(signif_TSS)>0 & signif_threshold[isim]>0) {
   models_TSS[[isim]]=toString(colnames(x)[signif_TSS])
  

} else {failTSS=failTSS+1;}

} else 

                   {regrdata3=as.data.frame(cbind(y,x[,signif]))
glm3=glm(y ~ ., data =regrdata3,family=binomial)
bestgl3=do.call("glmulti", list(glm3, method="h", level=1, crit="aic",plotty=F,report=F,family=family))
bestreg3=bestgl3@objects[[1]]
keepx3=names(bestreg3$model)[-1]
idsel3=match(keepx3, names(regrdata2))-1

signif_TSS=signif[idsel3]
  
}

} else {failTSS=failTSS+1;}



}


list(signif_TSS=signif_TSS,signif_threshold=signif_threshold,signif=signif,failTSS=failTSS)

}

prostate=tss(x=std(prostate[,-1]),y=prostate[,1],nrepeat=50, nrepeat5=10, nsim=1, percent_split=.5, cutoff=0.625, family="binomial")

toString(colnames(x)[prostate$signif_TSS])



