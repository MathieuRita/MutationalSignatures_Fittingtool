probabilities_to_exposure=function(mut_mat,
                                   signatures,
                                   crible,
                                   nboot_fit=1000)

    {

    
    final_expo=0*crible
    deviation_expo=0*crible
    
    
    # Fit each sample with the active signatures
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
    
    # Colnames and rownames of the exposure matrix
    if(length(colnames(mut_mat))==ncol(final_expo)){
        colnames(final_expo)=colnames(mut_mat)
    }
    
    if(length(colnames(signatures))==nrow(final_expo)){
        rownames(final_expo)=colnames(signatures)
    }
    
    return(list(final_expo,deviation_expo))

}




