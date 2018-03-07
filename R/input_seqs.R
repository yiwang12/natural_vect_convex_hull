library(seqinr) # read.fasta, s2c
# f_mat_dom convert all fragments of one seq into a matrix
f_mat_dom <-function(seq_inp,l_dom_inp,list_aa_inp){
aalist = list_aa_inp
Seq1 = seq_inp
l_dom = l_dom_inp
Seq1_1 = s2c(Seq1)
Seq1_num = unclass(factor(Seq1_1, levels = s2c(aalist)))
Seq1_num = as.numeric(Seq1_num)
n_dom = length(Seq1_num)-l_dom+1 # number of domain for each length
mat_d_all = array(dim=c(sum(n_dom),max(l_dom)))
list_length=array(dim=sum(n_dom)) # each sequence length

for(j in 1:length(l_dom)){
  mat_d = array(dim=c(n_dom[j],l_dom[j]))
  for(i in 1:n_dom[j]){
    mat_d[i,] = Seq1_num[i:(i+l_dom[j]-1)]}
  mat_d_all[ (sum(n_dom[1:j-1])+1) :sum(n_dom[1:j]) , 1:l_dom[j]] = mat_d
  list_length[ (sum(n_dom[1:j-1])+1) :sum(n_dom[1:j])]=l_dom[j]
}
  outpu = cbind(mat_d_all,list_length)
return(outpu)
}


