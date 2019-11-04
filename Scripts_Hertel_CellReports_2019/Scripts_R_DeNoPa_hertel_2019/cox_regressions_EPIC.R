#######################################
#### function to compute Cox-regression

cox_reg = function(dat, sdate, idate, outc, exp, conf){
  
  ## 'dat'   -- data set containing all data needed
  ## 'sdate' -- name of the variable of baseline examination
  ## 'idate' -- date of event
  ## 'outc'  -- name of the variable containing incident cases
  ## 'exp'   -- vector of exposure variables
  ## 'conf'  -- confounders as formula to be considered
  
  options(stringsAsFactors = F)
  require(survival)
  
  ## how many variables
  jj = length(exp)
  
  ## reorder the outc to begin with outcome with the highest number of valid observations
  exp = sapply(exp, function(x) sum(!is.na(dat[,x])))
  exp = names(sort(exp, decreasing = T))
  
  ## compute follow-up time
  dat$fol = as.numeric((as.POSIXct(dat[,idate], format = "%y-%m-%d") - as.POSIXct(dat[,sdate], format = "%y-%m-%d")))/365.25 
  
  print(summary(dat$fol))
  
  cat("run model:\n")
  
  ## loop over all exposure variables
  for(j in 1:jj){
    
    ## prepare formula for lm model
    ff = paste0("Surv(fol,",outc,")~ ", exp[j], " + ", conf)
    
    ## print to screen
    if(j == 1){
      cat(ff,"\n")
    }

    
    ## devide and prepare storage of results
    if(j == 1){
      
      ## cave: first outcome has to fit a valid model!!
      mm = coxph(as.formula(ff), data=dat, ties="breslow")
      ## get the summary
      ss = summary(mm)
      ## how many terms in the model
      dd = nrow(ss$coefficients)
      ## store names of coefficients
      nm = rownames(ss$coefficients)
      ## indicatior for saving
      ll = ss$nevent
      
      ## create storage of results
      beta           = array(data=NA, dim=c(jj, dd))
      colnames(beta) = paste("beta", nm, sep="_")
      se             = array(data=NA, dim=c(jj, dd))
      colnames(se)   = paste("se", nm, sep="_")
      pval           = array(data=NA, dim=c(jj, dd))
      colnames(pval) = paste("pval", nm, sep="_")
      # hr             = array(data=NA, dim=c(jj, dd))
      # colnames(hr)   = paste("hr", nm, sep="_")
      # cil            = array(data=NA, dim=c(jj, dd))
      # colnames(cil)  = paste("hr_cil", nm, sep="_")
      # ciu            = array(data=NA, dim=c(jj, dd))
      # colnames(ciu)  = paste("hr_ciu", nm, sep="_")
      n              = array(data=NA, dim=jj)
      nevent         = array(data=NA, dim=jj)
      
    }else{
      
      ## test if function will run
      ll = tryCatch({
        
        ## cave: first outcome has to fit a valid model!!
        mm = coxph(as.formula(ff), data=dat, ties="breslow")
        ## get the summary
        ss = summary(mm)
        ## indicatior for saving
        nn = ss$nevent
        
      }, error=function(e){
        cat("no compuation on ", exp[j], "\n")
        return(NA)
      }, warning=function(e){
        cat("no compuation on ", exp[j], "\n")
        return(NA)
      })
    }
    ## store if nothing failed
    if(!is.na(ll) & dd == nrow(ss$coefficients)){
      beta[j,]  = ss$coefficients[,1]
      # hr[j,]    = ss$coefficients[,2]
      pval[j,]  = ss$coefficients[,5]
      se[j,]    = ss$coefficients[,3]
      # cil[j,]   = ss$conf.int[,3]
      # ciu[j,]   = ss$conf.int[,4]
      n[j]      = ss$n
      nevent[j] = ss$nevent
    }
  }

  
  ## prepare data frame for return (suppress intercept)
  res = data.frame(exposure=exp, beta, se, pval, n, nevent)
  
  ## replace metabolite name with exposure in the names
  names(res) = gsub("M[0-9]*[0-9]", outc, names(res))

  
  return(res)
  
}