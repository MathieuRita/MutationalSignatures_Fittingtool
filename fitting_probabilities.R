fitting_probabilities=function(mut_mat,
                               signatures,
                               ind_prior_knowledge,
                               nboot=2,
                               initial_cutoff=0.2,
                               cutoff_remove=0.2,
                               cutoff_add=0.1,
                               cutoff_cs=0.01,
                               noise="Gaussian",
                               num_ite=2,
                               n_ite=10)
    
{
    
    # First initialization
    
    print("Vote 1")
    
    crible=sparse_filter(mut_mat=mut_mat,
                         signatures=signatures,
                         ind_prior_knowledge=ind_prior_knowledge,
                         nboot=nboot,
                         initial_cutoff=initial_cutoff,
                         cutoff_add=cutoff_add,
                         cutoff_remove=cutoff_remove,
                         cutoff_cs=cutoff_cs,
                         noise=noise,
                         num_ite=num_ite)
    
    for(i in seq(1,n_ite-1)){
        
        # Estimate the active signatures for several iterations
        print(paste("##### Vote",toString(i)))
        crible=crible+sparse_filter(mut_mat=mut_mat,
                                    signatures=signatures,
                                    ind_prior_knowledge=ind_prior_knowledge,
                                    nboot=nboot,
                                    initial_cutoff=initial_cutoff,
                                    cutoff_add=cutoff_add,
                                    cutoff_remove=cutoff_remove,
                                    cutoff_cs=cutoff_cs,
                                    noise=noise,
                                    num_ite=num_ite)
    }
    
    # For a signature, the activity probability is equal to the frequence of votes
    
    fitting_probabilities=crible/n_ite
    
    colnames(fitting_probabilities)=colnames(mut_mat)
    rownames(fitting_probabilities)=colnames(signatures)
    
    return(fitting_probabilities)
    
}

