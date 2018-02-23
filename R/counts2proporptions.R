# converts count data to proportions
#
# counts A data frame of pollen counts of size N_sites by N_taxa
# returns a data frame of the same size 
# Non-zero rows sum to 1
# Zero rows remain 0

counts2proportions <- function(counts){

  return(t(apply(counts, 1, function(x) { if (sum(x)>0) {x/sum(x)} else {0} })))
  
}