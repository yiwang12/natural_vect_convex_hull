library(quadprog)
library(Matrix)
source("function_dist_cal_domain_v2.R")
URL="/home/wangyi/revise/sum_frag_sep/list_nv_fragments_C1_cpkc2.RData"
load(URL)
load("!/inf_mat_all.RData")

dist_all_domain_fragments = list()

for(i in 1:516){
  dist_all_domain_fragments[[i]]=list()
  for(j in 1:417)
  {
    dist_all_domain_fragments[[i]][[j]]=list()
  }
}


#Ex. cPKC C1
load("~/nv_domain_C1_cpkc_mat.RData")

#> lenrage_C1_cpkc
#[1] 116 141
seq_range = c(1:516)

for( i in seq_range[which(!seq_range %in% inf_mat_all)]){
  dist_all_domain_fragments[[i]] = dist_cal_one_seq_all_fragments(nv_domain_C1_cpkc_mat,list_nv_fragments_C1_cpkc2[[i]],c(116:141))
}

dist_all_domain_C1_cpkc2 = dist_all_domain_fragments
save(dist_all_domain_C1_cpkc2,file="dist_all_domain_C1_cpkc2.RData")
