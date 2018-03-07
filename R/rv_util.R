library(moments)
# cal nvs for all fragments in ith seq
moment_nonzero <-function(vect_inp){
  vec2 = vect_inp[which(vect_inp!=0)]
  moment(vec2-mean(vec2)-1 ,order=2)
}

nvs_cal <-function(mat_d_all_inp){
mat_allseq = mat_d_all_inp[,1:ncol(mat_d_all_inp)-1]
list_length=mat_d_all_inp[,ncol(mat_d_all_inp)]
mat_allseq[is.na(mat_allseq)] = 0
n_seq=nrow(mat_allseq) # number of protein sequences
len_seq = ncol(mat_allseq)  # maximal length of sequences


vect=c(1:len_seq) # vect: 1 2 3 ..
mat = t(array(vect,dim=c(len_seq,n_seq)))
mat_all = array(mat, dim=c(dim(mat),20))

for (i in 1:20){
  mat_all[,,i][which(mat_allseq!=i)]=0
}


D1= apply(mat_all, c(1,3), moment_nonzero)

nt=list_length
#nt = array(dim=n_seq) #number of total aa
#for (s in 1:n_seq){
#  nt[s] = min(which(mat_allseq[s,]==0))-1
#}

nk=array(dim=c(n_seq,20)) # number of each type aa
mu=array(dim=c(n_seq,20)) # mean location of each aa
for (s in 1:n_seq){
  for ( i in 1:20){
    nk[s,i] = sum((mat_all[s,,i]!=0))
    mu[s,i]  = mean(which(mat_all[s,,i]!=0))
}}
mu[is.nan(mu)]=0

D2=D1
D=array(dim=dim(D2))

for  (i in 1:20){
D[,i]= D2[,i]/ nt}
D[is.nan(D)]=0
nv = cbind(nk,mu,D)
return(nv)
}

