remove_signatures=function(samples,
                           signatures,
                           set_of_signature,
                           cutoff=0.01)
  {
  
  "
  Takes the current set of signatures and remove the signatures that do not contribute to the reconstruction score.
  If removing a signature does not decrease the score more than (1-cutoff)*current_score, the signature is removed from the current set of signatures
  
  Args:
  - samples: matrix of samples you want to fit [nrow=96,ncol=num_of_samples]
  - signatures: matrix of all the signatures [nrow=96,ncol=num_of_signatures]
  - set_of_signature: vector with the indices of the current signatures
  - cutoff: cutoff to choose whether a signature is removed from the set of signatures (RSS criterion)
  
  Returns:
  - set_of_signature: the updated set of signatures
  "
  
  #BASELINE SCORE
  boot_1=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature]))
  baseline=RSS(rowMeans(samples),rowMeans(boot_1$reconstructed))
  
  new_best_score=0
  while(new_best_score<cutoff){
    
    #This condition is to prevent from a bug -> fitting with no signature (2 can be replaced by 1)
    if(length(set_of_signature)>2){
      
      #REMOVE STEPS
      RSSs=c()
      for(i in seq(1,length(set_of_signature))){
        boot_1=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature[-i]]))
        RSSs=c(RSSs,RSS(rowMeans(samples),rowMeans(boot_1$reconstructed)))
      }
      
      new_best_score=min((RSSs-baseline)/baseline)
      
      if(new_best_score<cutoff){
        set_of_signature=set_of_signature[-which.min(RSSs)]
      }
    }
    else{
      new_best_score=1
    }
  }
  
  return(set_of_signature)
}