# protein sequence classification using natural vector and convex hull method

### Transform fasta data to a list containing protein family, full sequence and domain sequence information
get_info_ncbi.R;	functions_get_info.R

### Using natural vector method, trasfer sequence data to 60-dimentional vector for each protein
input_seqs.R;nv_cal.R;nv_util.R

### Calculate distance from each protein (60-dimentional vector) to each protein family (convex hull)
dist_cal.R;function_dist_cal_domain.R
