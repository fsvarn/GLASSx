##################################################
# Perform focal LOH test for HLA as described in McGranahan et al 2017
# Author: Frederick Varn
# Date: 2021.11.29
# Reported in text, no figure
##################################################

load("/Users/varnf/Downloads/titan_dloh_probs.RData")

# This is your probability of loss, I just randomly simulated it here
probs <- noncodel_titan_dloh_probs

# Each iteration of this loop runs step 2 for a different sample
# The 10,000 permutations takes place in the loop by setting size = 10,000
# The perm_count variable has a slot for each iteration and sums the number of samples that are positive as the loop goes on

perm_count <- rep(0, 10000)
for(i in 1:length(probs))
{
	
	num_pos <- sample(x = c(0, 1), size = 10000, replace = TRUE, prob = c(1 - probs[i], probs[i]))
	perm_count <- perm_count + num_pos
}

# Now calculate a p-value
# Quantile = 0.95 for a one-sided test
one_sided_p <- quantile(perm_count,0.95)

# Quantile = 0.975 for a two-sided test
two_sided_p <- quantile(perm_count,0.975)

# Calculate your real p-value:
real_value <- 12		# Use this variable to set what the real proportion of samples with HLA LOH is
real_p <- (1 - (ecdf(perm_count)(real_value))) 

# Alternatively (same value): 
sum(perm_count > real_value)/ length(perm_count)



load("/Users/varnf/Downloads/titan_dloh_probs.RData")

# This is your probability of loss, I just randomly simulated it here
probs <- wt_titan_dloh_probs

# Each iteration of this loop runs step 2 for a different sample
# The 10,000 permutations takes place in the loop by setting size = 10,000
# The perm_count variable has a slot for each iteration and sums the number of samples that are positive as the loop goes on

perm_count <- rep(0, 10000)
for(i in 1:length(probs))
{
	
	num_pos <- sample(x = c(0, 1), size = 10000, replace = TRUE, prob = c(1 - probs[i], probs[i]))
	perm_count <- perm_count + num_pos
}

# Now calculate a p-value
# Quantile = 0.95 for a one-sided test
one_sided_p <- quantile(perm_count,0.95)

# Quantile = 0.975 for a two-sided test
two_sided_p <- quantile(perm_count,0.975)

# Calculate your real p-value:
real_value <- 21		# Use this variable to set what the real proportion of samples with HLA LOH is
real_p <- (1 - (ecdf(perm_count)(real_value))) 

# Alternatively (same value): 
sum(perm_count > real_value)/ length(perm_count)