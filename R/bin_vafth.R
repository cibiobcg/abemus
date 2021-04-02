bin_vafth <- function(bin,vaf,covbin,spec,replica=10){

  w <- vaf[which(covbin == bin)]

  vaf.thbin <- data.frame(bin=bin,
                          spec=spec,
                          th=as.numeric(quantile(w,probs = spec,na.rm = TRUE)),
                          n=length(w[!is.na(w)]),
                          n_vaf_gtz=length(w[w > 0]),
                          n_vaf_etz=length(w[w == 0]),
                          stringsAsFactors = FALSE)

  if(FALSE){
    vaf.thbin <- data.frame(bin=bin,
                            spec=spec,
                            th=as.numeric(quantile(w,probs = spec,na.rm = TRUE)),
                            keep=NA,
                            run=NA,
                            stringsAsFactors = FALSE)

    if( length(w) >= 100 ){

      nn <- length(w) - round(length(w) * seq(0.01,0.99,0.01))

      for(sz in nn){

        for(run in 1:replica){
          this <- data.frame(bin=bin,
                             spec=spec,
                             th=as.numeric(quantile(sample(w,size = sz, replace = FALSE),probs = spec,na.rm = TRUE)),
                             keep=sz,
                             run=run,
                             stringsAsFactors = FALSE)

          vaf.thbin <- rbind(vaf.thbin,this)
        }

      }

      vaf.thbin_full <- vaf.thbin %>%
        filter(is.na(keep)) %>%
        select(-keep,-run)

      vaf.thbin_subs <- vaf.thbin %>%
        filter(!is.na(keep)) %>%
        group_by(spec,keep) %>%
        summarise(cvar=sd(th,na.rm=TRUE)/mean(th,na.rm=TRUE),median=median(th,na.rm = TRUE)) %>%
        arrange(desc(keep))

      vaf.thbin <- full_join(x = vaf.thbin_full,y = vaf.thbin_subs,by='spec') %>%
        mutate(delta = abs(median - th))

      # tested only with one spec and one bin
      p <- ggplot(vaf.thbin, aes(keep, delta)) + geom_bar(stat = 'identity') +
        ggtitle(bin) +
        geom_hline(yintercept = 0.001,linetype="dashed", color = "red")
      print(p)

      return(vaf.thbin)
    }

  }

  return(vaf.thbin)
}
