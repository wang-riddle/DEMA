process_matrices <- function(matrix_list) {
  result_list <- list()
  for (i in seq_along(matrix_list)) {
    mat <- matrix_list[[i]]
    non_zero_cols <- which(colSums(mat != 0) > 0)
    filtered_mat <- mat[, non_zero_cols, drop = FALSE]
    positive_indices <- apply(filtered_mat, 2, function(col) {
      which(col > 0)
    })
    result_list[[i]] <- list(
      filtered_matrix = filtered_mat,
      positive_row_indices = positive_indices,
      original_col_indices = non_zero_cols 
    )
    
    cat(sprintf("matrix %d: primary=%d, after=%d\n", 
                i, ncol(mat), ncol(filtered_mat)))
  }
  
  return(result_list)
}

results=process_matrices(BETA)
VEC.list=list()
for (i in 1:length(BETA)) {
  MAT=BETA[[i]]
  MAT[MAT!=0]=1
  MAT[MAT==0]=0
  VEC=apply(MAT,1,sum)
  VEC.list[[i]]=VEC
}

## how many times the gene was selected in 100 times of trials
usage_import=matrix(0,8,p)
for (i in 1:8) {
  usage_import[i,]=VEC.list[[i]]
}
colnames(usage_import)=colnames(X1)
row.names(usage_import)=c("DEMA","LASSO","LASSOCOM",
                          "SCAD","SCADCOM","MCP",
                          "MCPCOM","TL")

non_zero_col=colSums(usage_import!=0)>0
usage_nonzero=usage_import[,non_zero_col,drop=F]
usage_nonzero_bin=ifelse(usage_nonzero>0,1,0)
##total gene numbers that are selected
apply(usage_nonzero_bin, 1, sum)

write.csv(usage_nonzero,"usage_nonzero_LUSC.csv")

########the order of gene times of each method
result_list_ord <- list()


for (i in 1:nrow(usage_import)) {
  sorted_indices <- order(usage_import[i, ], decreasing = TRUE)
  sorted_colnames <- colnames(usage_import)[sorted_indices]
  result_list_ord[[i]] <- sorted_colnames
}

names(result_list_ord) <- rownames(usage_import)

process_matrices <- function(matrix_list) {
  result_list <- list()
  
  for (i in seq_along(matrix_list)) {
    mat <- matrix_list[[i]]
    
    non_zero_cols <- which(colSums(mat != 0) > 0)
    filtered_mat <- mat[, non_zero_cols, drop = FALSE]
    positive_indices <- apply(filtered_mat, 2, function(col) {
      which(col > 0)
    })
    result_list[[i]] <- list(
      filtered_matrix = filtered_mat,
      positive_row_indices = positive_indices,
      original_col_indices = non_zero_cols  
    )
    cat(sprintf("matrix %d: primary=%d, after=%d\n", 
                i, ncol(mat), ncol(filtered_mat)))
  }
  
  return(result_list)
}

results=process_matrices(BETA)
VEC.list=list()
for (i in 1:length(BETA)) {
  MAT=BETA[[i]]
  MAT[MAT!=0]=1
  MAT[MAT==0]=0
  VEC=apply(MAT,1,sum)
  VEC.list[[i]]=VEC
}

## how many times the gene was selected in 100 times of trials
usage_import=matrix(0,8,p)
for (i in 1:8) {
  usage_import[i,]=VEC.list[[i]]
}
colnames(usage_import)=colnames(X1)
row.names(usage_import)=c("DEMA","LASSO","LASSOCOM",
                          "SCAD","SCADCOM","MCP",
                          "MCPCOM","TL")

non_zero_col=colSums(usage_import!=0)>0
usage_nonzero=usage_import[,non_zero_col,drop=F]
usage_nonzero_bin=ifelse(usage_nonzero>0,1,0)
##total gene numbers that are selected
apply(usage_nonzero_bin, 1, sum)

write.csv(usage_nonzero,"usage_nonzero_LUSC.csv")

########the order of gene times of each method
result_list_ord <- list()

for (i in 1:nrow(usage_import)) {
  sorted_indices <- order(usage_import[i, ], decreasing = TRUE)
  sorted_colnames <- colnames(usage_import)[sorted_indices]
  result_list_ord[[i]] <- sorted_colnames
}


names(result_list_ord) <- rownames(usage_import)

