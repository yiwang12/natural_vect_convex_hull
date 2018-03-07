get_location2 <- function(key_in,content_in){
  k=0
  out=array(dim=length(content_in))
  for (i in 1:length(content_in)){
    line = strsplit(content_gp_all[i], " ")
    if (length(line[[1]]   )>0){
    if (line[[1]][1]==key_in){k=k+1; out[k]=i}
  }}
  out=out[!is.na(out)]
  return(out)
}

get_location <- function(key_in,content_in){
  k=0
  out=array(dim=length(content_in))
  for (i in 1:length(content_in)){
    if (content_gp_all[i]==key_in){k=k+1; out[k]=i}
  }
  out=out[!is.na(out)]
  return(out)
}


get_location_grep<-function(key_in,content_in){
  k=0
  out=array(dim=length(content_in))
  for (i in 1:length(content_in)){
    t=content_in[i]
    a=grep(key_in,t)
    if(length(a)!=0)
    {k=k+1; out [k]=i}
  }
  out=out[!is.na(out)]
  return(out)
}


get_num <- function(loc_key,loc_num){
  out=array(dim=length(loc_key))
  for (i in 1:length(loc_key)){
    a=min(which(loc_num>loc_key[i] ))  # which seq num the "C2" belong to
    out[i]=a
  }
  return(out)
}


get_fam_sameline <- function(loc_key,loc_num){
  out=array(dim=length(loc_key))
  for (i in 1:length(loc_key)){
    a=which(loc_num==loc_key[i]) # which seq num the "C2" belong to
    out[i]=a
  }
  return(unique(out))
}




get_acc_list <- function(loc_acc_in,content_in){
  out = list()
  for (i in 1:length(loc_acc_in)){
    a=content_in[loc_acc_in[i]]
    b=strsplit(a, "\\s+")[[1]]
    out[[i]] = b[2]
  }
  return(out)
}

get_SOURCE_list <- function(loc_SOURCE_in,content_in){
  out = list()
  for (i in 1:length(loc_SOURCE_in)){
    a=content_in[loc_SOURCE_in[i]]
    b=strsplit(a, "     ")[[1]]
    out[[i]] = b[2]
  }
  return(out)
}

extract_domain_loc<-function(loc_domain_in,loc_region_in,content_in){
  out=list()
  out_st = array(dim=length(loc_domain_in))
  out_en = array(dim=length(loc_domain_in))
  for (i in 1:length(loc_domain_in)){
    a = loc_domain_in[i]
    b=loc_region_in[max(which(loc_region_in<a))]
    c=content_in[b]
    d=strsplit(c, "\\s+")[[1]]
    e=d[3]
    info=strsplit(e,"..",fixed=TRUE)[1]
    out_st[i] = as.numeric(info[[1]][1])
    out_en[i] = as.numeric(info[[1]][2])
  }
  out$st = out_st
  out$en = out_en
  return(out)
}

getdomain_Seq <- function(list_acc,pkc_fa,loc_domain_st,loc_domain_en,num_seqf_inp){
  outpu = list()
  for (i in 1:length(list_acc)){
    num_f = num_seqf_inp[i]
    seq= as.character(pkc_fasta_all[num_f])
    seq2 = s2c(seq)
    loc_domain_ = c(loc_domain_st[i],loc_domain_en[i])
    if(!is.na(loc_domain_[1]) &  !is.na(loc_domain_[2] )){
      seq3 = seq2[loc_domain_[1]:loc_domain_[2]]
      seq4=paste(seq3, sep="", collapse="")
      outpu[[i]] = seq4
    }
  }
  return(outpu)
}



read_file_content<-function(URL_file_in){
  file_in <- file(URL_file_in,open="r")
  out <-readLines(file_in) #loc of "//" 
  close(file_in)
  return(out)
}

num_famname_xpkc_all_in=num_famname_npkc_all
get_table_acc_sour<-function(num_famname_xpkc_all_in,acc_list_all_in,list_SOURCE_in){
  a=list()
  for(i in 1:length(num_famname_xpkc_all_in)){
    a[[i]]=list()
    seq_num = num_famname_xpkc_all_in[i]
    a[[i]][[1]] = seq_num
    a[[i]][[2]] = as.character(list_SOURCE_in[seq_num] [[1]] [1])
  }
  b=matrix(unlist(a), nrow = length(num_famname_xpkc_all_in), byrow = TRUE)
  return(b)
}


