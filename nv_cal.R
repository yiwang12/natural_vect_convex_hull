library(seqinr) # read.fasta, s2c
library(quadprog) # cal dist
library(moments) # cal nv
source("rv_util.R") # cal nv !set
source("input_seqs.R") # cal nv !set

pkc_samp <- read.fasta(file = "/home/wangyi/polytopes/Human_kinase_protein_516-2.fasta", as.string = TRUE, seqtype = "AA")
aa_list = "GAVLIFWYDHNEKQMRSTCP"


l_domain = c(125:140) #length range of kinase domain 

shi_nvs = list()

for(i in 51:100){
  if(nchar(pkc_samp[i])>max(l_domain)){
    seq= as.character(pkc_samp[i])
    mati = f_mat_dom(seq,l_domain,aa_list) # convert all fragments of one seq into a matrix
    nvs = nvs_cal(mati)  # cal nvs for all fragments in ith seq
    shi_nvs[[i]]=nvs
    save(shi_nvs,file="shi_nvs.RData")}
}



