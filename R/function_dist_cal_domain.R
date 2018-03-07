

library(quadprog)
library(Matrix)

dist_onegroup<-function(X_n,Y){
  X_n <- X_n[!duplicated(X_n), ]
  n_vec = dim(X_n)[1]
  H = X_n %*% t(X_n)
 # H2=H
  H <- nearPD(H)$mat
  f =  Y %*% t(X_n)  
  # equalities
  A.eq = t(array(1,dim=n_vec))
  b.eq = 1
  # inequalities
  #non
  # lower-bounds 
  A.lbs = diag(1,n_vec,n_vec)
  b.lbs = array(0,dim=n_vec)
  # upper-bounds 
  A.ubs = -diag(1,n_vec,n_vec)
  b.ubs = -array(1,dim=n_vec)
  sol <- solve.QP(Dmat = H,
                  dvec = f,
                  Amat = t(rbind(A.eq, A.lbs, A.ubs)),
                  #Amat = t(rbind(A.eq, A.ge, A.lbs, A.ubs)),
                  bvec = c(b.eq,  b.lbs, b.ubs),
                  #bvec = c(b.eq, b.ge, b.lbs, b.ubs),
                  meq = 1,factorized=FALSE)
  
  dist= sqrt((sol$value + sum(Y^2)/2)*2) 
  return(dist)
}



dist_cal_one_seq_all_fragments<-function(X_n_inpu,Ys,domain_lenrange_i){
  dist_mat = list()
    for(i in domain_lenrange_i){
      {
        dist_mat[[i]] = list()
        for(j in 1:length(Ys[[i]])){
          dist_mat[[i]][[j]] = dist_onegroup(X_n_inpu,Ys[[i]][[j]])
        }
      }
  }
  return(dist_mat)
}






