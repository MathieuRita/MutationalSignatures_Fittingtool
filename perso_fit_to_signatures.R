perso_fit_to_signatures=function(mut_mat,
                                 signatures,
                                 ind_prior_knowledge,
                                 n_ite=10,
                                 prob_cutoff=0.3,
                                 nboot_train=2,
                                 nboot_fit=100,
                                 initial_cutoff=0.2,
                                 cutoff=0.01,
                                 cutoff_add=0.1){
    
    sparse_filter = function(mut_mat,
    signatures,
    ind_prior_knowledge,
    nboot=2,
    initial_cutoff=0.2,
    cutoff_add=0.1,
    cutoff_cs=0.01,
    noise="Gaussian",
    num_ite=2)
    
    #Select the set of signatures by sample using several layers
    crible=sparse_filter(mut_mat,signatures,ind_prior_knowledge,nboot=nboot_train,initial_cutoff,cutoff,cutoff_add)
    
    for(i in seq(1,n_ite-1)){
        print(print(paste("##### Vote",toString(i))))
        crible=crible+sparse_filter(mut_mat,signatures,ind_prior_knowledge,nboot=nboot_train,initial_cutoff,cutoff,cutoff_add)
    }
    
    return(crible/n_ite)
   

fitting_probabilities=function(mut_mat,
                               signatures,
                               ind_prior_knowledge,
                               nboot=2,
                               initial_cutoff=0.2,
                               cutoff_add=0.1,
                               cutoff_cs=0.01,
                               noise="Gaussian",
                               num_ite=2)
    
{
    
    # First initialization
    crible=sparse_filter(mut_mat=mut_mat,
                         signatures=signatures,
                         ind_prior_knowledge=ind_prior_knowledge,
                         nboot=nboot,
                         initial_cutoff=initial_cutoff,
                         cutoff_add=cutoff_add,
                         cutoff_cs=cutoff_cs,
                         noise=noise,
                         num_ite=num_ite)
    
    for(i in seq(1,n_ite-1)){
        
        # Estimate the active signatures for several iterations
        print(print(paste("##### Vote",toString(i))))
        crible=crible+sparse_filter(mut_mat=mut_mat,
                                    signatures=signatures,
                                    ind_prior_knowledge=ind_prior_knowledge,
                                    nboot=nboot,
                                    initial_cutoff=initial_cutoff,
                                    cutoff_add=cutoff_add,
                                    cutoff_cs=cutoff_cs,
                                    noise=noise,
                                    num_ite=num_ite)
    }
    
    # For a signature, the activity probability is equal to the frequence of votes
    
    fitting_probabilities=crible/n_ite
    
    return(fitting_probabilities)
    
}



   
'
    
    
# TO DO

probabilities_to_exposure=function(mut_mat,
                                   signatures,
                                   crible,
)

    
    final_expo=0*crible
    deviation_expo=0*crible
    
    
    
    
    # Bootstrap version
    for(i in seq(1,ncol(mut_mat))){
        ind_sg=which(crible[,i]==1)
        fboot=fit_boot(mut_mat[,i],signatures[,ind_sg],sample_name=colnames(mut_mat[i]),nboot=nboot_fit)
        
        sample_expo=rowMeans(fboot$contribution)
        
        sample_std=rep(0,nrow(fboot$contribution))
        for(k in seq(1,nrow(fboot$contribution))){
            sample_std[k]=std.error(fboot$contribution[k,])
        }
        
        for(j in seq(1,length(sample_expo))){
            final_expo[ind_sg[j],i]=sample_expo[j]
            deviation_expo[ind_sg[j],i]=sample_std[j]
        }
    }
    
    if(length(colnames(mut_mat))==ncol(final_expo)){
        colnames(final_expo)=colnames(mut_mat)
    }
    
    if(length(colnames(signatures))==nrow(final_expo)){
        rownames(final_expo)=colnames(signatures)
    }
    
    
    return(list(final_expo,deviation_expo))
    '
}




