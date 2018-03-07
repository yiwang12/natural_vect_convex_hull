source("functions_get_info.R")
## sum loc of each seq
URL_file = "~/File_combined.txt"
content_gp_all = read_file_content(URL_file)
loc_num_all=get_location('//',content_gp_all)

#> length(loc_num_all)
#[1] 3047
##### sum loc of domain
loc_c1_all = get_location_grep('/region_name="C1',content_gp_all)
#length(loc_c1_all)
#[1] 3053
loc_c2_all = get_location_grep('/region_name="C2',content_gp_all)
#> length(loc_c2_all)
#[1] 929
loc_kin_all = get_location_grep('/region_name="STKc_',content_gp_all)
#> length(loc_kin_all)
#[1] 1327
loc_pb1_all = get_location_grep('/region_name="PB1',content_gp_all)
#length(loc_pb1_all)
#[1] 332

### get correponding seq number of each key word 
len_c2= length(loc_c2_all)
len_c1= length(loc_c1_all)
len_k= length(loc_kin_all)
len_p= length(loc_pb1_all)


li_c2_num_all = array(dim=len_c2)
li_c1_num_all = array(dim=len_c1)
li_kin_num_all = array(dim=len_k)
li_p_num_all = array(dim=len_p)


li_c2_num_all = get_num(loc_c2_all,loc_num_all)

#> dim(li_c2_num_all)
#[1] 929
li_c1_num_all = get_num(loc_c1_all,loc_num_all)

#> dim(li_c1_num_all)
#[1] 3053
li_kin_num_all = get_num(loc_kin_all,loc_num_all)

#> dim(li_kin_num_all)
#[1] 1327

li_p_num_all = get_num(loc_pb1_all,loc_num_all)

#> length(li_p_num_all)
#[1] 332


#### accesion loc
loc_acc_all =get_location_grep('ACCESSION',content_gp_all) 

#> length(loc_acc_all)
#[1] 3047

len_acc_all = length(loc_acc_all)

li_acc_num_all = get_num(loc_acc_all,loc_num_all)

#> sum(li_acc_num_all!=c(1:len_acc_all))
#[1] 0

acc_list_all = get_acc_list(loc_acc_all,content_gp_all)

#> length(acc_list_all)
#[1] 3047


####extract domain  sequence location
loc_region_all = get_location_grep('Region',content_gp_all)
  
#> length(loc_region_all)
#[1] 8157

#extract C2 sequence location
C2_domain_loc = extract_domain_loc(loc_c2_all,loc_region_all,content_gp_all)
C2_St_all= C2_domain_loc$st
C2_En_all =  C2_domain_loc$en
#!!some loc info not available NA exist


#extract kinase sequence location

k_domain_loc = extract_domain_loc(loc_kin_all,loc_region_all,content_gp_all)
k_St_all= k_domain_loc$st
k_En_all =  k_domain_loc$en

#extract pb1 sequence location

pb1_domain_loc = extract_domain_loc(loc_pb1_all,loc_region_all,content_gp_all)
pb1_St_all= pb1_domain_loc$st
pb1_En_all =  pb1_domain_loc$en
##! some missing info, NAs introduced

#extract c1 sequence location
c1_domain_loc = extract_domain_loc(loc_c1_all,loc_region_all,content_gp_all)
C1_St_all= c1_domain_loc$st
C1_En_all =  c1_domain_loc$en

##! some missing info, NAs introduced



## merge all info
li_info_all = list()

for(i in 1:length(acc_list_all)){ #c2
  num = which(li_c2_num_all==i)
  if(length(num)>0){
    li_info_all$c2_st_all[i] =  C2_St_all[num] 
    li_info_all$c2_en_all[i] =  C2_En_all[num]
  }
}
##! some missing info, NAs introduced


#> length(which(!is.na(li_info_all$c2_en_all)))
#[1] 923
#> length(which(!is.na(li_info_all$c2_st_all)))
#[1] 907


for(i in 1:length(acc_list_all)){ #pb1
  num = which(li_p_num_all==i)
  if(length(num)>0){
    li_info_all$pb1_st_all[i] =  pb1_St_all[num] 
    li_info_all$pb1_en_all[i] =  pb1_En_all[num]
  }
}

#> length(which(!is.na(li_info_all$pb1_st_all)))
#[1] 309
#> length(which(!is.na(li_info_all$pb1_en_all)))
#[1] 326

for(i in 1:length(acc_list_all)){ #kinase
  num = which(li_kin_num_all==i)
  if(length(num)>0){
    li_info_all$k_st_all[i] =  k_St_all[num] 
    li_info_all$k_en_all[i] =  k_En_all[num]
  }
}

#> length(which(!is.na(li_info_all$k_st_all)))
#[1] 1327
#> length(which(!is.na(li_info_all$k_en_all)))
#[1] 1327



for(i in 1:length(acc_list_all)){ #c1
  num = which(li_c1_num_all==i)
  if(length(num)==2){
    li_info_all$c1_st_all[i] =  C1_St_all[num[1]] 
    li_info_all$c1_en_all[i] =  C1_En_all[num[2]]
  }
}

#> length(which(!is.na(li_info_all$c1_st_all)))
#[1] 1227
#> length(which(!is.na(li_info_all$c1_en_all)))
#[1] 1235

for(i in 1:length(acc_list_all)){ #c1
  num = which(li_c1_num_all==i)
  if(length(num)==1){
    li_info_all$c1_st_all_apkc[i] =  C1_St_all[num]
    li_info_all$c1_en_all_apkc[i] =  C1_En_all[num]
  }
}

#> length(which(!is.na(li_info_all$c1_st_all)))
#[1] 1227
#> length(which(!is.na(li_info_all$c1_en_all)))
#[1] 1235



### match with fasta through accesion num
URL_fasta = "/Users/yiwang/Dropbox/Mac/submit_revise_v3/NCBI_data/NCBI_sum_ac/whseq/sequence.fasta.txt"
content_fasta = read_file_content(URL_fasta)

loc_seqf_all = get_location_grep('>',content_fasta)
  
num_seqf_all = array(dim=length(acc_list_all))# the number in fasta file for the ith accesion 
for (i in 1:length(acc_list_all)){
  #t=file_fa[i]
  a=grep(acc_list_all[i],content_fasta)
  if(length(a)==1){
    num_seqf_all[i] = which(loc_seqf_all==a) # running slow, 3min
  }
}

i=981

#save(num_seqf_all,file="num_seqf_all.RData")

#### scan for domain seq from fasta file using matched number of accesion
library(seqinr) # read.fasta, s2c
URL_fasta = "/Users/yiwang/Dropbox/Mac/submit_revise_v3/NCBI_data/NCBI_sum_ac/whseq/sequence.fasta.txt"
pkc_fasta_all= read.fasta(file = URL_fasta, as.string = TRUE, seqtype = "AA")


#num_seqf_inp=num_seqf_all
#list_acc=acc_list_all
#pkc_fa=pkc_fasta_all
#loc_domain_st=li_info_all$k_st_all
#loc_domain_en=li_info_all$k_en_all


domain_k_all=getdomain_Seq(acc_list_all,pkc_fasta_all, li_info_all$k_st_all, li_info_all$k_en_all,num_seqf_all)#NULL exist
domain_c1_all=getdomain_Seq(acc_list_all,pkc_fasta_all, li_info_all$c1_st_all, li_info_all$c1_en_all,num_seqf_all)#NULL exist
domain_c2_all=getdomain_Seq(acc_list_all,pkc_fasta_all, li_info_all$c2_st_all, li_info_all$c2_en_all,num_seqf_all) #NULL exist
domain_pb1_all=getdomain_Seq(acc_list_all,pkc_fasta_all, li_info_all$pb1_st_all, li_info_all$pb1_en_all,num_seqf_all) #NULL exist


#### PKC family classify
PKC_family = c("conventional|Conventional|cPKC|alpha|beta|gamma","novel|Novel|nPKC|theta|epsilon|delta| eta","atypical|Atypical|aPKC|lambda|iota|zeta")

loc_famname_cpkc_all = get_location_grep(PKC_family[1],content_fasta)

loc_famname_npkc_all = get_location_grep(PKC_family[2],content_fasta)

loc_famname_apkc_all = get_location_grep(PKC_family[3],content_fasta)


#num the class loc

num_famname_cpkc_all = get_fam_sameline(loc_famname_cpkc_all,loc_seqf_all)
num_famname_npkc_all = get_fam_sameline(loc_famname_npkc_all,loc_seqf_all) 
num_famname_apkc_all =get_fam_sameline(loc_famname_apkc_all,loc_seqf_all)



#> length(num_famname_cpkc_all)
#[1] 531
#> length(num_famname_npkc_all)
#[1] 887
#> length(num_famname_apkc_all)
#[1] 687

#> intersect(num_famname_cpkc_all,num_famname_npkc_all)
#integer(0)
#> intersect(num_famname_cpkc_all,num_famname_apkc_all)
#[1] 1628 !! remove this seq !!
#> intersect(num_famname_npkc_all,num_famname_apkc_all)
#integer(0)

### remove seqs with intercept between different families
1628
num_accession_intercept = which(num_seqf_all==1628) #num_seqf_all : the number in fasta file for the ith accesion 
#> num_accession_intercept
#[1] 527

### remove non-pkc seqs by checking name 
PKC_remove = "like|interacting|partial|related|containing|inhibitor|similar|region|similarity"

loc_famname_cpkc_all = get_location_grep(PKC_remove,content_fasta)

num_accession_name = which(num_seqf_all %in% loc_famname_cpkc_all) #num_seqf_all : the number in fasta file for the ith accesion 


### remove seqs without reference
loc_ref_all = get_location_grep("REFERENCE",content_gp_all)
num_ref_all = array(dim=length(loc_ref_all))
for ( i in 1:length(loc_ref_all)){
  num_ref_all[i]= min(which(loc_num_all>loc_ref_all[i]))
}
num_ref_2_all = unique(num_ref_all) #seqs number with reference

num_accession_remove_ref_intercept = num_ref_2_all[which(num_ref_2_all!=527)]

num_accession_remove_ref_intercept_name =num_accession_remove_ref_intercept[!num_accession_remove_ref_intercept %in% num_accession_name]

length(num_accession_remove_ref_intercept)
#[1] 1822


length(num_accession_remove_ref_intercept_name)
#[1] 1819



num_seqf_all2 = unique(num_seqf_all)
num_seqf_all_filtered = num_seqf_all2[which(num_seqf_all2 %in% num_accession_remove_ref_intercept_name)]

#> length(num_seqf_all2)
#[1] 3047
#> length(num_seqf_all)
#[1] 3047

pkc_filtered = pkc_fasta_all[num_seqf_all_filtered]


num_seqf_all_filtered
length(num_seqf_all_filtered)


#> length(pkc_filtered)
#[1] 1818


#> sum(is.na(num_seqf_all))
#[1] 0

#############domain corresponding to pkc_filtered

num_famname_cpkc_all_accesion_order_filtered = which( num_seqf_all  %in% num_seqf_all_filtered[which (num_seqf_all_filtered %in% num_famname_cpkc_all)])
domain_k_cpkc_filtered = domain_k_all[num_famname_cpkc_all_accesion_order_filtered] 
domain_C1_cpkc_filtered = domain_c1_all[num_famname_cpkc_all_accesion_order_filtered] 
domain_C2_cpkc_filtered = domain_c2_all[num_famname_cpkc_all_accesion_order_filtered] 

#length(num_famname_cpkc_all_accesion_order_filtered)
#[1] 324
#> sum((domain_k_cpkc_filtered)!="NULL")
#[1] 214
#> sum((domain_C1_cpkc_filtered)!="NULL")
#[1] 218
#> sum((domain_C2_cpkc_filtered)!="NULL")
#[1] 281

num_famname_npkc_all_accesion_order_filtered =  which( num_seqf_all  %in% num_seqf_all_filtered[which (num_seqf_all_filtered %in% num_famname_npkc_all)])
domain_k_npkc_filtered = domain_k_all[num_famname_npkc_all_accesion_order_filtered] 
domain_C1_npkc_filtered = domain_c1_all[num_famname_npkc_all_accesion_order_filtered] 
domain_C2_npkc_filtered = domain_c2_all[num_famname_npkc_all_accesion_order_filtered] 

#length(num_famname_npkc_all_accesion_order_filtered)
#[1] 538
#> sum((domain_k_npkc_filtered)!="NULL")
#[1] 247
#> sum((domain_C1_npkc_filtered)!="NULL")
#[1] 423
#> sum((domain_C2_npkc_filtered)!="NULL")
#[1] 125

num_famname_apkc_all_accesion_order_filtered =  which( num_seqf_all  %in% num_seqf_all_filtered[which (num_seqf_all_filtered %in% num_famname_apkc_all)])
domain_k_apkc_filtered = domain_k_all[num_famname_apkc_all_accesion_order_filtered] 
domain_C1_apkc_filtered = domain_c1_all_apkc[num_famname_apkc_all_accesion_order_filtered] 
#domain_C2_apkc_filtered = domain_c2_all[num_famname_apkc_all_accesion_order_filtered] 
domain_C2_pb1_filtered = domain_pb1_all[num_famname_apkc_all_accesion_order_filtered] 


#length(num_famname_apkc_all_accesion_order_filtered)
#[1] 440
#> sum((domain_k_apkc_filtered)!="NULL")
#[1] 227
#> sum((domain_C1_apkc_filtered)!="NULL")
#[1] 326
#> sum((domain_C2_pb1_filtered)!="NULL")
#[1] 209

#> sum((domain_C2_apkc_filtered)!="NULL")
#[1] 1
#> which((domain_C2_apkc_filtered)!="NULL")
#[1] 117 
#> num_famname_apkc_all_accesion_order_filtered[117]
#[1] 422
#> num_seqf_all[422]
#[1] 623
 
domain_C1_cpkc_filtered
domain_k_all=getdomain_Seq(acc_list_all,pkc_fasta_all, li_info_all$k_st_all, li_info_all$k_en_all,num_seqf_all)#NULL exist
domain_c1_all=getdomain_Seq(acc_list_all,pkc_fasta_all, li_info_all$c1_st_all, li_info_all$c1_en_all,num_seqf_all)#NULL exist
domain_c2_all=getdomain_Seq(acc_list_all,pkc_fasta_all, li_info_all$c2_st_all, li_info_all$c2_en_all,num_seqf_all) #NULL exist
domain_pb1_all=getdomain_Seq(acc_list_all,pkc_fasta_all, li_info_all$pb1_st_all, li_info_all$pb1_en_all,num_seqf_all) #NULL exist
domain_c1_all_apkc=getdomain_Seq(acc_list_all,pkc_fasta_all, li_info_all$c1_st_all_apkc, li_info_all$c1_en_all_apkc,num_seqf_all) #NULL exist


num_famname_cpkc_all = get_fam_sameline(loc_famname_cpkc_all,loc_seqf_all)
num_famname_npkc_all = get_fam_sameline(loc_famname_npkc_all,loc_seqf_all) 
num_famname_apkc_all =get_fam_sameline(loc_famname_apkc_all,loc_seqf_all)


