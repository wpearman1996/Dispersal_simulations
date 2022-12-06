forward_mod <- function (initial, prob = 0, d = 1, gens = 150, keep = FALSE, 
    pool = NULL, limit.sim = FALSE, coeff.lim.sim = 1, sigm = 0.1, 
    filt = NULL, prob.death = NULL, method.dist = "euclidean", 
    plot_gens = FALSE, iteration = NULL, sub_iteration = NULL, first.gen = 0, stage = NULL) 
{
    if (!is.numeric(prob) | prob < 0) {
        stop("Probability of migration or mutation must be a number belonging to ", 
            "[0; 1] interval.")
    }
    if (!is.numeric(d) | d < 0) {
        stop("Number of individuals that die in each time step must be a positive", 
            " number.")
    }
    if (!is.numeric(gens) | gens <= 0) {
        stop("Number of generations must be a positive number.")
    }
    if (!is.logical(keep)) {
        stop("keep parameter must be a boolean.")
    }
    if (!is.logical(limit.sim)) {
        stop("limiting similarity parameter must be a boolean.")
    }
    if (!is.numeric(coeff.lim.sim)) {
        stop("coeff.lim.sim parameter must be numeric.")
    }
    if (!is.numeric(sigm) | sigm < 0) {
        stop("sigm parameter must be a positive number.")
    }
    if ((method.dist %in% c("euclidean", "maximum", "manhattan", 
        "canberra", "binary", "minkowski")) == FALSE) {
        stop("Provided distance does not exist. See stats::dist function for help.")
    }
    if (!is.logical(plot_gens)) {
        stop("plot_gens parameter must be a boolean.")
    }
    if ((is.character(initial) | is.vector(initial)) & (limit.sim | 
        !is.null(filt))) {
        stop("Trait information must be provided along with species identity in", 
            " the initial community for niche - based dynamics")
    }
    if (!is.matrix(initial) & !is.data.frame(initial) & (limit.sim | 
        !is.null(filt))) {
        stop("Misdefined initial community")
    }
    if (!limit.sim & is.null(filt)) {
        message("Simulation of a neutral community")
    }
    if (is.character(pool)) {
        pool <- data.frame(id = 1:length(pool), sp = pool, trait = rep(NA, 
            length(pool)), stringsAsFactors = FALSE)
        if (limit.sim | !is.null(filt)) {
            message("No trait information provided in the regional pool")
            pool[, 3] <- runif(nrow(pool))
            message("Random trait values attributed to individuals of the regional", 
                " pool")
            colnames(pool) <- c("id", "sp", "trait")
        }
    }
    if (!is.null(pool)) {
        if (ncol(pool) < 2) {
            stop("The regional pool is misdefined (at least two columns ", 
                "required when a matrix or data frame is provided)")
        }
        else if (ncol(pool) == 2) {
            message("No trait information provided in the regional pool")
        }
        if ((!limit.sim | !is.null(filt)) & ncol(pool) < 3) {
            pool[, 3] <- runif(nrow(pool))
            message("Random (uniform) trait values attributed to individuals of ", 
                "the regional pool")
        }
        colnames(pool) <- c("id", "sp", "trait")
    }
    if (is.character(initial)) {
        J <- length(initial)
        init_comm <- data.frame(id = paste("init", 1:J, sep = ""), 
            sp = initial, trait = rep(NA, J), stringsAsFactors = FALSE)
    }
    else {
        if (ncol(initial) < 3) {
            message("Two-column initial community: assumed to represent species ", 
                "and trait information; individual ids will be generated")
            J <- nrow(initial)
            init_comm <- data.frame(id = paste("init", 1:J, sep = ""), 
                sp = initial[, 1], trait = initial[, 2], stringsAsFactors = FALSE)
        }
        else {
            J <- nrow(initial)
            init_comm <- initial
        }
    }
    if (J < d) 
        stop("The number of dead individuals per time step ", 
            "cannot be greater than community size")
    if ((limit.sim | !is.null(filt)) & any(is.na(init_comm[, 
        3]))) {
        stop("Trait information must be provided in initial community ", 
            "composition for niche-based dynamics")
    }
    colnames(init_comm) <- c("id", "sp", "trait")
    new.index <- 0
    next_comm <- init_comm
    sp_t <- c()
    ind_t <- c()
    dist.t <- c()
    if (keep) {
        comm_through_time <- c() #empty placeholder
        
        ifelse(is.null(iteration), dirnam <- paste0("output/iteration_", Sys.Date()), dirnam <- paste0("output/iteration_", formatC(iteration, width = 5, flag = 0)))
        
        if(!is.null(sub_iteration)){dirnam <- paste0(dirnam, "/", formatC(sub_iteration, width = 3, flag = 0))}
        dir.create(dirnam, showWarnings=F) #instead of making a list in memory output to textfiles - first create a folder
    }
    for (i in 1:gens) {
        if (keep) {
        	if(i == gens){comm_through_time <- next_comm}
            if(is.null(sub_iteration)){write.csv(next_comm, paste0(dirnam, "/generation_", formatC(i, width = 5, flag = "0", format = "d"), "_", stage, ".csv"))}else{
            	gen_mod <- first.gen + i
            	write.csv(next_comm, paste0(dirnam, "/generation_", formatC(gen_mod, width = 5, flag = "0", format = "d"), "_", stage, ".csv"))
            } #save each step as named csv
        }
        next_comm <- pick(next_comm, d = d, prob = prob, pool = pool, 
            prob.death = prob.death, limit.sim = limit.sim, coeff.lim.sim = coeff.lim.sim, 
            sigm = sigm, filt = ifelse(is.list(filt), filt[[i]], filt), new.index = new.index, 
            method.dist = "euclidean")
        sp_t <- c(sp_t, length(unique(next_comm$com$sp)))
        ind_t <- c(ind_t, length(unique(next_comm$com$ind)))
        if (!is.null(limit.sim)) {
            dist.t <- c(dist.t, next_comm$dist.t)
        }
        new.index <- next_comm$new.index
        next_comm <- next_comm$com
    }
    if (plot_gens) {
        uniq_df <- data.frame(gens = 1:gens, ind_t = ind_t, sp_t = sp_t, 
            stringsAsFactors = FALSE)
        if (requireNamespace("ggplot2", quietly = TRUE)) {
            plot_individuals <- ggplot(uniq_df, aes_string("gens", 
                "ind_t")) + geom_line() + geom_line(aes_string("gens", 
                "ind_t"), size = 1) + labs(x = "Number of generations", 
                y = "Number of distinct ancestors")
            plot_species <- ggplot(uniq_df, aes_string("gens", 
                "sp_t")) + geom_line(size = 1) + labs(x = "Number of generations", 
                y = "Number of species")
            print(plot_individuals)
            print(plot_species)
        }
    }
    if (!is.null(limit.sim)) {
        if (keep) 
            return(list(com_t = comm_through_time, sp_t = sp_t, 
                dist.t = dist.t, pool = pool))
        else return(list(com = next_comm, sp_t = sp_t, dist.t = dist.t, 
            pool = pool))
    }
    else {
        if (keep) 
            return(list(com_t = comm_through_time, sp_t = sp_t, 
                pool = pool))
        else return(list(com = next_comm, sp_t = sp_t, pool = pool))
    }
}
