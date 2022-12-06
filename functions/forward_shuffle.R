forward_shuffle <- function (init_pools = NULL, prob = 0, d = 1, gens = 150, keep = TRUE, 
    pool = NULL, limit.sim = FALSE, coeff.lim.sim = 1, sigm = 0.1, 
    filt = NULL, prob.death = NULL, method.dist = "euclidean", 
    plot_gens = FALSE, cohort.size = 10, shuffle.freq = 10, shuffle.prop = 0.1, iteration = NULL, first.gen = 0, stage = NULL){ 
    	#init_pools = list of initial pools, one for each cohort - if NULL create pools using coalescence
    	#cohort.size = number of co-migrants in cohort 
    	#shuffle.freq = shuffle every x generations 
    	#shuffle.prop = [0-1] proportion of biome from each to shuffle
		#iteration = used for naming folders		
		
     #outputcall
     ifelse(is.null(iteration), dirnam <- paste0("output/iteration_", Sys.Date()), dirnam <- paste0("output/iteration_", formatC(iteration, width = 5, flag = 0)))
     dir.create(dirnam, showWarnings=F) #output to textfiles - first create a folder
     filt_info <- ifelse(is.list(filt), {paste(filt_dispersal_seq[1], filt_dispersal_seq[length(filt)/2], filt_dispersal_seq[length(filt)], sep = " -> ")}, {get(strsplit(strsplit(deparse(filt)[2], "\\(")[[1]][2], "\\,")[[1]][1])})

     suppressWarnings(write.table(data.frame(STAGE = stage[1], iteration = iteration[1], prob = prob[1], gens = gens[1], filter = filt_info, limit.sim = limit.sim[1], coeff.lim.sim = coeff.lim.sim[1], sigm = sigm[1], prob.death = prob.death[1]), paste0(dirnam, "/_call.csv"), row.names = FALSE, append = TRUE, sep = ","))
     #outputcall
	
	
	if(!any(c(is.null(cohort.size), is.null(shuffle.freq), is.null(shuffle.prop)))){
	
	n_shuffles <- ceiling(gens/shuffle.freq) 
	
	if(is.null(init_pools)){
			init_pools <- lapply(1:cohort.size, function(x){coalesc(J = 10000, m = prob, filt = filt, pool = pool)$com})
								}					
		
	#shuffle
	for(l in 1:n_shuffles){
		if(is.list(filt)){
			filt_subset <- lapply((((l-1)*shuffle.freq)+1):(l* shuffle.freq), function(x){filt[[x]]}) #get list of niche filters for this subset from overarching dispersal filter list
		}else{filt_subset <- filt}

		if(l == 1){current_pools <- lapply(1:length(init_pools), function(initx){
				forward_mod(initial = init_pools[[initx]], prob = prob, prob.death = prob.death, gens = shuffle.freq, pool = pool, 
         	    filt = filt_subset, limit.sim = limit.sim, coeff.lim.sim = coeff.lim.sim, sigm = sigm, keep = keep, iteration = iteration, sub_iteration = initx, first.gen = first.gen, stage = stage)$com
				})}else{
				
		ind2shuffle <- sample(1:nrow(current_pools[[1]]), size = (nrow(current_pools[[1]]) * shuffle.prop))
		shuffle_pool <- do.call(rbind, lapply(current_pools, function(x){
			x[ind2shuffle,]
		}))
		
		shuffle_pool$sp <- sapply(shuffle_pool$sp, function(x){
			ifelse(substr(x, nchar(x)-2, nchar(x)) == "_sh", return(x), return(paste0(x, "_sh"))) #labelling shuffled/shared species
				})
				
		
		current_pools <- lapply(current_pools, function(x){
			pool <- x[-ind2shuffle,]
			pool <- rbind(pool, shuffle_pool[sample(nrow(shuffle_pool), size = (nrow(current_pools[[1]]) * shuffle.prop)), ])
			return(pool)		
		})
		
		current_pools <- lapply(1:length(current_pools), function(initx){
						forward_mod(initial = current_pools[[initx]], prob = prob, prob.death = prob.death, gens = shuffle.freq, pool = pool, 
                    	       filt = filt_subset, limit.sim = limit.sim, coeff.lim.sim = coeff.lim.sim, sigm = sigm, keep = keep, iteration = iteration, sub_iteration = initx, first.gen = first.gen+((l-1)*shuffle.freq), stage = stage)$com				
					})}
					}}else{
				current_pools <- lapply(1:length(init_pools), function(initx){
				filt_subset <- filt	
				forward_mod(initial = init_pools[[initx]], prob = prob, prob.death = prob.death, gens = gens, pool = pool, 
         	    filt = filt_subset, limit.sim = limit.sim, coeff.lim.sim = coeff.lim.sim, sigm = sigm, keep = keep, iteration = iteration, sub_iteration = initx, first.gen = first.gen, stage = stage)$com
					})}
		return(current_pools)			
					}