############################################
#### function to compute Cox-regression meta
#### analyses for EPIC

meta_cox_reg = function(dat, sdate, idate, outc, exp, conf, label){
  
  ## 'dat'   -- data set containing all data needed
  ## 'sdate' -- name of the variable of baseline examination
  ## 'idate' -- date of event
  ## 'outc'  -- name of the variable containing incident cases
  ## 'exp'   -- vector of exposure variables
  ## 'conf'  -- confounders to be considered
  ## 'label' -- data frame to ensure metabolite identification
  
  ## package for meta-analysis
  require(metafor)
  
  ## function for regression
  source("../scripts/cox_reg.R")
  
  ## first batch
  res.1 = cox_reg(subset(dat, BatchMTBL == 2), sdate, idate, outc, exp, conf)
  ## skip NA
  res.1 = na.omit(res.1)
  cat("found ", nrow(res.1), " associations in the first batch\n")

  ## second batch
  res.2 = cox_reg(subset(dat, BatchMTBL == 3), sdate, idate, outc, exp, conf)
  ## skip NA
  res.2 = na.omit(res.2)
  cat("found ", nrow(res.2), " associations in the second batch\n")

  ## merge results
  res   = merge(res.1[, c("exposure", "n", "nevent", grep(outc, names(res.1), value=T))],
                res.2[, c("exposure", "n", "nevent", grep(outc, names(res.2), value=T))],
                by=c("exposure"), all=T)
  
  ## how many observations in common
  res$n      = res$n.x + res$n.y
  res$nevent = res$nevent.x + res$nevent.y
  
  ## at least 10 events at total
  res = subset(res, nevent > 10)
  cat("found ", nrow(res), " associations at all\n")

  
  ## define columns for meta analyses
  # vars = c(paste0("beta_", outc, c(".x",".y")), paste0("se_", outc, c(".x",".y")))
  # 
  # ## run meta analyses
  # res[, c("beta.meta", "se.meta", "pval.meta", "heterogenity", "pval.hetero")] = t(apply(res[, vars], 1, 
  #                                                                                        function(x){a = rma(yi=x[1:2], sei = x[3:4], method="FE");
  #                                                                                        return(c(a$beta, a$se, a$pval, a$I2, a$QEp))}))
  
  ## define how many different betas were found
  
  beta <- grep(paste("beta.*", outc, sep="_"), names(res.1), value=T)
  print(beta)
  
  if(length(beta) == 1){
    ## get ride of the proxy
    beta <- gsub("beta_", "", beta)
    ## define columns for meta analyses
    vars = c(paste0("beta_", beta, c(".x",".y")), paste0("se_", beta, c(".x",".y")))
    
    ## run meta analyses
    res[, c("beta.meta", "se.meta", "pval.meta", "heterogenity", "pval.hetero")] = t(apply(res[, vars], 1, 
                                                                                           function(x){a = rma(yi=x[1:2], sei = x[3:4], method="FE");
                                                                                           return(c(a$beta, a$se, a$pval, a$I2, a$QEp))}))
    
  }else{
    for(k in gsub("beta_", "", beta)){
      ## define columns for meta analyses
      vars = c(paste0("beta_", k, c(".x",".y")), paste0("se_", k, c(".x",".y")))
      
      ## run meta analyses
      res[, paste(c("beta.meta", "se.meta", "pval.meta", "heterogenity", "pval.hetero"), k, sep=".")] = t(apply(res[, vars], 1, 
                                                                                                                function(x){a = rma(yi=x[1:2], sei = x[3:4], method="FE");
                                                                                                                return(c(a$beta, a$se, a$pval, a$I2, a$QEp))}))
    } 
  }
  
  
  
  
  ## add metabolite identification
  res = merge(res, label, by.x="exposure", by.y="ID", all.x=T)
  
  return(res)
}