remove_specific_signatures=function(samples,
                                    signatures,
                                    set_of_signature,
                                    cutoff=0.01)
  {
  
  
  
  #BASELINE SCORE
  boot_1=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature]))
  baseline=RSS(rowMeans(samples),rowMeans(boot_1$reconstructed))
  
  new_best_score=0
  
  continue=T
  while(continue){
    #This condition is to prevent from a bug -> fitting with no signature (2 can be replaced by 1)
    continue=F
    if(length(set_of_signature)>2){
      #REMOVE STEPS
      RSSs=c()
      stds=c()
      for(i in seq(1,length(set_of_signature))){
        boot_1=fit_to_signatures(samples,as.matrix(signatures[,set_of_signature[-i]]))
        RSSs=c(RSSs,RSS(rowMeans(samples),rowMeans(boot_1$reconstructed)))
        stds=c(stds,std.error(signatures[,set_of_signature[i]]))
      }
      
      ind=which.min(RSSs*stds)
      
      new_best_score=(RSSs[ind]-baseline)/baseline
      
      
      if(new_best_score<cutoff){
        set_of_signature=set_of_signature[-ind]
        continue=T
      }
    }
  }
  
  return(set_of_signature)
}