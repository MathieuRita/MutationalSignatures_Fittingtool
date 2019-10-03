sparse_filter = function(mut_mat,
                         signatures,
                         ind_prior_knowledge,
                         nboot=2,
                         initial_cutoff=0.2,
                         cutoff_add=0.1,
                         cutoff_cs=0.01,
                         cutoff_remove=0.2,
                         noise="Gaussian",
                         num_ite=2)

    {
        
    "
    Takes the current set of signatures and remove the signatures that do not contribute to the reconstruction score.
    If removing a signature does not decrease the score more than (1-cutoff)*current_score, the signature is removed from the current set of signatures
        
    Args:
    - mut_mat: matrix of samples you want to fit [nrow=96,ncol=num_of_samples]
    - signatures: matrix of all the signatures [nrow=96,ncol=num_of_signatures]
    - ind_prior_knowledge: vector of the indices of signatures we expect to be present in the set of samples
    - initial_cutoff: cutoff for the initial layer
    - cutoff: cutoff to choose whether a signature is removed from the set of signatures (RSS criterion)
    - cutoff_add: cutoff for the adding layer
    - cutoff_remove: cutoff for the removing layer
    - cutoff_exposure: cutoff for the final layer
    - noise: Gaussian or Poisson: for the resampling, select the type of noise added to the data
        
    Returns:
    - sieve: for each samples, a binary list (0 or 1) indicating which signature is active in a sample
    "
    
    # Initialization of the sieve
    sieve=matrix(0L,ncol=ncol(mut_mat),nrow=ncol(signatures))
    
    # For each sample, find the set of active signatures
    
    for(i in seq(1,ncol(mut_mat))){
        
        print(i)

        # Step 0 : Resampling (dependingo on the type of noise the user choose)

        if(noise=="Gaussian"){
            
            sample=mut_mat[,i]
            samples=matrix(0L,ncol=nboot,nrow=96)
            samples[,1]=sample
            for(j in seq(2,nboot)){
                samples[,j]=profile_bootstrap(sample)
            }
            
            }
        
        if(noise=="Poisson"){
            
            sample=mut_mat[,i]
            samples=matrix(0L,ncol=nboot,nrow=96)
            for(j in seq(1,nboot)){
                samples[,j]=add_poisson_noise(sample)
            }
        }
        
        # Step 1 : Initial set of signatures
        # To draw the initial set of signatures, we select the most active signature when we fit the sample with all the signatures. We concatenate them with the set of expected signatures (prior_knowledge)
        
        
        fit_add=fit_to_signatures(samples,as.matrix(signatures))
        set_of_sig_to_add=as.numeric(which(rowMeans(fit_add$contribution)>initial_cutoff*sum(sample)))
        set_of_signature=union(ind_prior_knowledge,set_of_sig_to_add)
        
        # Step 2 : Perform the adding and removing steps iteratively
        
        for(j in seq(1,num_ite)){
            #Removing layer
            #final_set_of_signature=remove_specific_signatures(samples,signatures,set_of_signature,cutoff=200*cutoff)
            set_of_signature=remove_signatures(samples,signatures,set_of_signature,cutoff=cutoff_remove)
            
            # Adding layer
            set_of_signature=add_signatures(samples,signatures,set_of_signature,cutoff_add=0.5*cutoff_add,cutoff_cs=0.2*cutoff_cs)
        }

        # Step 3: Eventually remove the signatures with small activities
        set_of_signature=remove_small_exposure(samples,signatures,set_of_signature,cutoff_exposure=0.05)
        
        # Step 4: Write in the sieve the list of active signatures
        sieve[set_of_signature,i]=1
        
    }
    
    return(sieve)
}
