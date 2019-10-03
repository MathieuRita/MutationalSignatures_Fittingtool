remove_small_exposure=function(samples,
                               signatures,
                               set_of_signature,
                               cutoff_exposure=0.05)
  
  {

  "
  As a final step, the function remove all the signatures that contribute less than a certain percentage of the total contribution
  
  Args:
  - samples: matrix of samples you want to fit [nrow=96,ncol=num_of_samples]
  - signatures: matrix of all the signatures [nrow=96,ncol=num_of_signatures]
  - set_of_signature: vector with the indices of the current signatures
  - cutoff_exposure: cutoff to choose whether a signature is removed from the set of signatures (RSS criterion)
  
  Returns:
  - set_of_signature: the updated set of signatures
  "
  
  #Total number of mutations in the sample
  total_num_mut=mean(colSums(samples))
  
  # loop_bool is a boolean that checks if we have to continue removing signatures
  loop_bool=T
  
  while(loop_bool){
    loop_bool=F
    
    if(length(set_of_signature)>1){
      
      # We fit or samples
      fit=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature]))
      
      #Compute the mean exposure by signature
      meanExpo=rowMeans(fit$contribution)
      
      # We remove the small activities one by one starting with the signature with the smallest activity
      min=min(meanExpo)
      
      # if the smallest activity accounts for less than cutoff_exposure*the total number of mutations in the sample, we remove it
      
      if(min<total_num_mut*cutoff_exposure){
        set_of_signature=set_of_signature[-which.min(meanExpo)]
        loop_bool=T
      } 
      
    }
    
  }
  
  return(set_of_signature)
  
}
