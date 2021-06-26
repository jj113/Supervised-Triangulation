library(rARPACK);library(MFPCA);library(tidyverse);library(MASS);require(zipfR);require(fOptions);require(Matrix)
library(pec);library(fda);library(survAUC);library(Triangulation);library(BPST)

source("utility_fun.R")

theta = 0.01 # an example for a pre-specified theta (needs to be tuned)
npc = 2 # number of basis 
ncr = 40 # dimension of the simulation images
lambda = c(0, 1, 10, 10^2, 10^3, 10^6) # tunning parameter for bernstein

train_dat.id = read_rds('train_dat')
test_dat.id = read_rds('test_dat')
tdat.id = train_dat.id

V.est = read_rds('V.est') # vertices for triangles for simulation setting
Tr.est = read_rds('Tr.est') # triangles
Z = read_rds('Z') # expanded grid points
d.est = 2; r = 1
ind.inside=inVT(V.est, Tr.est, Z[,1], Z[,2])$ind.inside

#----------- for training ---------------
sup_fun = sfpca_img(type = "bernstein", Y = train_dat.id$Z1, train_dat.id = train_dat.id,
                    theta = theta, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                    r = r, Z = Z, ncr = ncr)
  
o_tfs = as.matrix(sup_fun[[1]])

b.basis = as.matrix(sup_fun[[3]]); b.scores = as.matrix(sup_fun[[4]])

score_sup = t(o_tfs) %*% (b.basis%*%t(b.scores))

na.y = rep(NA, 40*40)
na.y[ind.inside] = o_tfs[,1]
image(matrix(na.y, nrow = 40, ncol = 40))

if(npc == 2){
  score_sup = t(score_sup)
}else{
  score_sup = cbind(score_sup[1,])
}
  
score_names = c()
for(q in 1:(ncol(score_sup))){
    tname = paste("score", as.character(q), sep = "")
    score_names = c(score_names, tname)
}

tdat.id = tdat.id[,c(1:4)]
tdat.id = cbind(tdat.id, score_sup)
  
colnames(tdat.id)[(ncol(tdat.id) - ncol(score_sup) + 1) : ncol(tdat.id)] = score_names
  
fmla = as.formula(paste("Surv(time,event) ~ ", paste(score_names, collapse= "+")))
  
fitted_obj = coxph(fmla, data = tdat.id, x = T, y = T) 

#------------------------- TEST ---------------------------------------
  
est = bernstein(Y = test_dat.id$Z1, V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                lambda = lambda)

b.basis.test = as.matrix(est[[1]]); b.scores.test = as.matrix(t(est[[2]]))

score_sup_test = t(o_tfs) %*% ((b.basis.test%*%t(b.scores.test)))

if(npc == 2){
    score_sup_test = t(score_sup_test)
}else{
    score_sup_test = cbind(score_sup_test[1,])
}
  
  
test_dat.id = test_dat.id[,c(1:4)]

test_dat.id = cbind(test_dat.id, score_sup_test)
  
colnames(test_dat.id)[(ncol(test_dat.id) - ncol(score_sup_test) + 1) : ncol(test_dat.id)] = score_names
  
train.prob = predict(fitted_obj)
surv.prob = predict(fitted_obj, newdata = test_dat.id)
  
Surv.train = Surv(tdat.id$time, tdat.id$event)
Surv.test = Surv(test_dat.id$time, test_dat.id$event)
  
times = seq(0, 15, by = 0.1)

#---- prediction accuracy -----#
# cumulative AUC 
auc.sFPCA = AUC.uno(Surv.train, Surv.test, surv.prob, times)$iauc
auc.sFPCA  

# cumulative BS
SB = pec::pec(object = fitted_obj, formula = fmla,
              traindata = tdat.id,
              data=test_dat.id)
  
SB.sFPCA = crps(SB)[1]
SB.sFPCA 
 
 
