matchedwake <- function(data, t_window, spat_window, treatment, control, dependent, matchColumns, t_unit = "days", estimation = "lm", weighted = FALSE, estimationControls = c(), TCM = FALSE, deleteSUTVA = FALSE, alpha1 = 0.05, alpha2 = 0.1, memory = 1, match.default = TRUE, ...){
  missing_arguments <- c()
  terminate <- FALSE
  if(missing(data)){   
    missing_arguments <- append(missing_arguments,'\n  data')
    terminate <- TRUE
  }
  if (!is.element('timestamp',names(data))){
    missing_arguments <- append(missing_arguments,'\n  data: timestamp column is missing')
    terminate <- TRUE
  }
  if (!is.element('lat',names(data))){
    missing_arguments <- append(missing_arguments,'\n  data: lat column is missing')
    terminate <- TRUE
  }
  if (!is.element('long',names(data))){
    missing_arguments <- append(missing_arguments,'\n  data: long column is missing')
    terminate <- TRUE
  }
  if (missing(t_window)){
    missing_arguments <- append(missing_arguments,'\n  t_window')
    terminate <- TRUE
  }
  if (missing(spat_window)){
    missing_arguments <- append(missing_arguments,'\n  spat_window')
    terminate <- TRUE
  }
  if (missing(treatment)){
    missing_arguments <- append(missing_arguments,'\n  treatment')
    terminate <- TRUE
  }
  if (missing(control)){
    missing_arguments <- append(missing_arguments,'\n  control')
    terminate <- TRUE
  }
  if (missing(dependent)){
    missing_arguments <- append(missing_arguments,'\n  dependent')
    terminate <- TRUE
  }
  if (missing(matchColumns)){
    missing_arguments <- append(missing_arguments,'\n  matchColumns')
    terminate <- TRUE
  }
  if (!is.data.frame(data)){
    missing_arguments <- append(missing_arguments,'\n  data')
    terminate <- TRUE
  }
  datecheck <- try(as.Date(data$timestamp, format= "%Y-%m-%d %H:%M:%S"))
  if(any(class(datecheck) == "try-error" || is.na(datecheck))){
    missing_arguments <- append(missing_arguments,'\n  timestamp (date format must be "YYYY-MM-DD hh:mm:ss)')
    terminate <- TRUE
  }
  if (!is.numeric(t_window) || length(t_window)!=3){
    missing_arguments <- append(missing_arguments,'\n  t_window')
    terminate <- TRUE
  }
  if (t_window[1]<2){
    missing_arguments <- append(missing_arguments,'\n  t_window (Note: minimum time window size is 2 in order to calculate a trend)')
    terminate <- TRUE
  }
  if (!(abs((t_window[2]-t_window[1])/t_window[3]) - round((t_window[2]-t_window[1])/t_window[3])  < .Machine$double.eps^0.5)){
    missing_arguments <- append(missing_arguments,'\n  t_window (non-integer division of interval)')
    terminate <- TRUE
  }
  if (!is.numeric(spat_window) || length(spat_window)!=3){
    missing_arguments <- append(missing_arguments,'\n  spat_window')
    terminate <- TRUE
  }
  if (!(abs((spat_window[2]-spat_window[1])/spat_window[3]) - round((spat_window[2]-spat_window[1])/spat_window[3])  < .Machine$double.eps^0.5)){
    missing_arguments <- append(missing_arguments,'\n  spat_window (non-integer division of interval)')
    terminate <- TRUE
  }
  if (!is.character(treatment[1]) || length(treatment)!=2){
    missing_arguments <- append(missing_arguments,'\n  treatment')
    terminate <- TRUE
  }
  if (!is.character(control[1]) || length(control)!=2){
    missing_arguments <- append(missing_arguments,'\n  control')
    terminate <- TRUE 
  }
  if (!is.character(dependent) || length(dependent)!=2){
    missing_arguments <- append(missing_arguments,'\n  dependent')
    terminate <- TRUE
  }
  if (!is.character(matchColumns) || length(matchColumns)<1){
    missing_arguments <- append(missing_arguments,'\n  matchColumns (Note: at least one matching covariate is required)')
    terminate <- TRUE
  }
  missing_arguments <- append(missing_arguments,'\n')
  if (terminate){
    message("Error in matchedwake(data, t_window, spat_window,...):")
    message("The following required arguments of matchedwake() are MISSING or MIS-SPECIFIED:", missing_arguments)  
    stop("matchedwake(data, t_window, spat_window,...) was not executed!", call.=FALSE)
  }else{
    call <- match.call()
    cat('MWA: Analyzing causal effects in spatiotemporal event data.\n\nPlease always cite:\nSchutte, S., Donnay, K. (2014). "Matched wake analysis: Finding causal relationships in spatiotemporal event data." Political Geography 41:1-10.\n\n\nATTENTION: Depending on the size of the dataset and the number of spatial and temporal windows, iterations can take some time!\n\n')
    cat(paste('matchedwake(data = ', call$data,', t_window = ', call$t_window,', spat_window = ', call$spat_window,', ...):\n', sep = "")) 
    wakes <- slidingWake(data, t_unit, t_window, spat_window, treatment, control, dependent, matchColumns, estimationControls, memory)
    if (deleteSUTVA){
      wakes <- subset(wakes,(wakes$SO_pre==0 & wakes$MO_pre==0))
    }
    cat('Estimating causal effect on matched data...')
    matchedwake <- slideWakeMatch(wakes, alpha1, matchColumns, estimation, weighted, estimationControls, TCM, match.default, ...)
    cat('                               [OK]\n')
    if (!match.default){
      cat("ATTENTION: match.default set to FALSE, data was not matched!\n")
    }
    matchedwake$parameters <- list(t_window = t_window, spat_window = spat_window, treatment = treatment, control = control, dependent = dependent, matchColumns = matchColumns, t_unit = t_unit, estimation = estimation, estimationControls = estimationControls, TCM = TCM, deleteSUTVA = deleteSUTVA, alpha1 = alpha1, alpha2 = alpha2, memory = memory, match.default = match.default, ...)
    matchedwake$call <- call    
    class(matchedwake) <- "matchedwake"
    cat('Analysis complete!\n\nUse summary() for an overview and plot() to illustrate the results graphically.\n')
    return(matchedwake)
  } 
}

summary.matchedwake <- function(object, detailed = FALSE, ...){
  significant_results <- subset(object$estimates,object$estimates$pvalue<=object$parameters$alpha1)
  significant_results$estimate <- round(significant_results$estimate,digits=3)
  significant_results$pvalue <- round(significant_results$pvalue,digits=3)
  row.names(significant_results) <- NULL
  if (object$parameters$estimation == "lm"){
    significant_results$adj.r.squared <- round(significant_results$adj.r.squared,digits=4)
    names(significant_results) <- c(paste("Time[",object$parameters$t_unit,"]",sep=""), "Space[km]","EffectSize", "p_value","adj.Rsquared")
  }else{
    names(significant_results) <- c(paste("Time[",object$parameters$t_unit,"]",sep=""), "Space[km]","EffectSize", "p_value")
  }
  if (detailed){
    significant_matching <- subset(object$matching,object$estimates$pvalue<=object$parameters$alpha1)
    significant_SUTVA <- subset(object$SUTVA,object$estimates$pvalue<=object$parameters$alpha1)
    significant_results$perc_treatment <- round(100*significant_matching$treatment_post/(significant_matching$control_post+significant_matching$treatment_post),digits=1)
    significant_results$L1_post <- round(significant_matching$L1_post,digits=3)
    significant_results$CommonSupport_post <- significant_matching$commonSupport_post
    significant_results$SO <- 100*significant_SUTVA$SO
    significant_results$MO <- 100*significant_SUTVA$MO
    if (object$parameters$estimation == "lm"){
      names(significant_results) <- c(paste("Time[",object$parameters$t_unit,"]",sep=""), "Space[km]", "EffectSize", "p_value", "adj.Rsquared", "%Treatment", "L1metric", "%supp", "%S0", "%MO")
    }else{
      names(significant_results) <- c(paste("Time[",object$parameters$t_unit,"]",sep=""), "Space[km]", "EffectSize", "p_value", "%Treatment", "L1metric", "%supp", "%S0", "%MO")
    }
  }
  cat("Results:\n")
  if (nrow(significant_results)>0){
    cat("The method has identified combinations of temporal and spatial window sizes with significant estimates. The table below gives the estimated effect sizes and p-values for these combinations.\n\n")
    if (detailed){
      cat("NOTE: if a detailed summary was requested, the table also provides summary post-matching statistics as well as statistics on SUTVA violations:\nMatching:\n - '%Treatment' is short for the percentage of treatment events; the closer this is to 50%, the better the balance of treatment and control cases after matching\n - 'L1metric' measures the covariate balance after matching\n - '%supp' is short for percentage of common support after matching, a standard measure for the overlap of the range of covariate values\nSUTVA violations:\n - '%SO' gives the percentage of cases in which two or more treatment (or control) episodes overlap\n - '%MO' provides the fraction of overlapping treatment and control episodes\n\n")
    }
    return(significant_results)
  }
  if (nrow(significant_results)==0){ 
    cat("The method could not identify any space and time windows with significant estimates.\n")
  }
}

plot.matchedwake <- function(x, ...){
  density <- 3
  lty <- 2
  lwd <- 2
  pdata <- x$estimates
  t_unit <- x$parameters$t_unit
  alpha1 <- x$parameters$alpha1
  alpha2 <- x$parameters$alpha2
  yAxis <- sort(unique(pdata$t_window))
  xAxis <- sort(unique(pdata$spat_window))
  pswcplot <- matrix(nrow=length(xAxis)+2,ncol=length(yAxis)+2)
  eswcplot <- matrix(nrow=length(xAxis)+2,ncol=length(yAxis)+2)
  for (y in 1:length(yAxis)){
    for (x in 1:length(xAxis)){
      eswcplot [x+1,y+1] <- as.numeric(pdata$estimate[pdata$t_window == yAxis[y] & pdata$spat_window == xAxis[x]])
      pswcplot [x+1,y+1] <- as.numeric(pdata$pvalue[pdata$t_window == yAxis[y] & pdata$spat_window == xAxis[x]])
    }
  }
  for (y in 1:(length(yAxis)+2)){
    eswcplot [1,y] <- eswcplot [2,y]
    eswcplot [length(xAxis)+2,y] <- eswcplot [length(xAxis)+1,y]
  }
  for (x in 1:(length(xAxis)+2)){
    eswcplot [x,1] <- eswcplot [x,2]
    eswcplot [x,length(yAxis)+2] <- eswcplot [x,length(yAxis)+1]
  }
  filled.contour(x = seq(-1/(length(xAxis)-1),1+1/(length(xAxis)-1),length.out=length(xAxis)+2),y = seq(-1/(length(yAxis)-1),1+1/(length(yAxis)-1),length.out=length(yAxis)+2),eswcplot,xlab = "Spatial window in kilometers",
                 ylab=paste("Temporal window in ",t_unit,sep=""),
                 nlevels = 20,
                 color.palette = gray.colors,
                 xlim = c(-1/(length(xAxis)-1)/2,1+1/(length(xAxis)-1)/2),
                 ylim = c(-1/(length(yAxis)-1)/2,1+1/(length(yAxis)-1)/2),
                 plot.axes= {
                   axis(1,labels=xAxis,at=seq(0,1,by=1/(length(xAxis)-1)))
                   axis(2,at=seq(0,1,by=1/(length(yAxis)-1)),labels=yAxis)
                   yhalfRectSize <- (1 / (length(yAxis) -1)) / 2
                   xhalfRectSize <- (1 / (length(xAxis) -1)) / 2
                   ystepsize <- 1 / (length(yAxis) - 1)
                   xstepsize <- 1 / (length(xAxis) - 1)
                   for (y in 1:length(yAxis)){
                     for(x in 1:length(xAxis)){
                       if (is.na(pswcplot[x+1,y+1]) | is.nan(pswcplot[x+1,y+1]) | pswcplot[x+1,y+1] > alpha2){
                         rect(c(((x - 1)*xstepsize) - xhalfRectSize,((x - 1)*xstepsize) - xhalfRectSize),c(((y - 1)*ystepsize) - yhalfRectSize,((y - 1)*ystepsize) - yhalfRectSize),c(((x - 1)*xstepsize) + xhalfRectSize,((x - 1)*xstepsize) + xhalfRectSize),c(((y - 1)*ystepsize) + yhalfRectSize,((y - 1)*ystepsize)+yhalfRectSize),density=density,lwd=lwd)
                       } else if (is.na(pswcplot[x+1,y+1]) | is.nan(pswcplot[x+1,y+1]) | pswcplot[x+1,y+1] > alpha1 & pswcplot[x+1,y+1] <= alpha2){
                         rect(c(((x - 1)*xstepsize) - xhalfRectSize,((x - 1)*xstepsize) - xhalfRectSize),c(((y - 1)*ystepsize) - yhalfRectSize,((y - 1)*ystepsize) - yhalfRectSize),c(((x - 1)*xstepsize) + xhalfRectSize,((x - 1)*xstepsize) + xhalfRectSize),c(((y - 1)*ystepsize) + yhalfRectSize,((y - 1)*ystepsize)+yhalfRectSize),density=density,lty=lty, lwd=lwd)
                       }
                     }
                   }
                 }
  )
}

print.matchedwake <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  print(summary(x))
}

slideWakeMatch <- function(wakes, alpha1, matchColumns, estimation, weighted, estimationControls, TCM, match.default, ...){
  t_window <- spat_window <- NULL 
  ceminputs <- list(...)
  terms <- which(c("cem.baseline.group","cem.datalist","cem.cutpoints","cem.grouping","cem.eval.imbalance","cem.k2k","cem.method","cem.mpower","cem.L1.breaks","cem.L1.grouping","cem.verbose","cem.keep.all") %in% names(ceminputs)) 
  ceminputs[which(!(names(ceminputs) %in% c("cem.baseline.group","cem.datalist","cem.cutpoints","cem.grouping","cem.eval.imbalance","cem.k2k","cem.method","cem.mpower","cem.L1.breaks","cem.L1.grouping","cem.verbose","cem.keep.all"))) ] <- NULL
  if (length(ceminputs) > 0){
    names(ceminputs) <- unlist(lapply(1:length(ceminputs), function(x) unlist(strsplit(names(ceminputs[x]),split="cem."))[2]))
  }
  lminputs <- list(...)
  terms <- which(c("lm.subset","lm.weights","lm.na.action","lm.method","lm.model","lm.x","lm.y","lm.qr","lm.singular.ok","lm.contrasts","lm.offset") %in% names(lminputs)) 
  lminputs[which(!(names(lminputs) %in% c("lm.subset","lm.weights","lm.na.action","lm.method","lm.model","lm.x","lm.y","lm.qr","lm.singular.ok","lm.contrasts","lm.offset")))] <- NULL
  if (length(lminputs) > 0){
    names(lminputs) <- unlist(lapply(1:length(lminputs), function(x) unlist(strsplit(names(lminputs[x]),split="lm."))[2]))
  }
  attinputs <- list(...)
  terms <- which(c("att.model","att.extrapolate","att.ntree","att.x","att.vars","att.object","att.plot","att.ecolors") %in% names(attinputs)) 
  attinputs[which(!(names(attinputs) %in% c("att.model","att.extrapolate","att.ntree","att.x","att.vars","att.object","att.plot","att.ecolors")))] <- NULL
  if (length(attinputs) > 0){
    names(attinputs) <- unlist(lapply(1:length(attinputs), function(x) unlist(strsplit(names(attinputs[x]),split="att."))[2]))
  }
  nbinputs <- list(...)
  terms <- which(c("glm.nb.weights","glm.nb.subset","glm.nb.na.action","glm.nb.start","glm.nb.etastart","glm.nb.mustart","glm.nb.control","glm.nb.model","glm.nb.x","glm.nb.y","glm.nb.contrasts","glm.nb.init.theta") %in% names(nbinputs)) 
  nbinputs[which(!(names(nbinputs) %in% c("glm.nb.weights","glm.nb.subset","glm.nb.na.action","glm.nb.start","glm.nb.etastart","glm.nb.mustart","glm.nb.control","glm.nb.model","glm.nb.x","glm.nb.y","glm.nb.contrasts","glm.nb.init.theta")))] <- NULL
  if (length(nbinputs) > 0){
    names(nbinputs) <- unlist(lapply(1:length(nbinputs), function(x) unlist(strsplit(names(nbinputs[x]),split="glm.nb."))[2]))
  }
  formula <- "dependent_post ~ dependent_pre + treatment"
  if (length(estimationControls) > 0){
    form_unlist <- unlist(strsplit(formula,split="\\+"))
    formula <- form_unlist[1]
    if (length(form_unlist) > 2){
      for (term in 2:(length(form_unlist)-1)){
        formula <- paste(formula," +",form_unlist[term])
      }
    }
    for (term in 1:length(estimationControls)){
      formula <- paste(formula,"+ ",estimationControls[term],sep="")
    }
    formula <- paste(formula," +",form_unlist[length(form_unlist)])
  }
  nterms <- length(unlist(strsplit(formula,split="\\+")))
  t_windows <- sort(unique(wakes$t_window))
  spat_windows <- sort(unique(wakes$spat_window))
  if (estimation == "lm"){
    estimates<- as.data.frame(matrix(nrow = (length(t_windows) * length(spat_windows)),ncol = 5))
    names(estimates) <- c("t_window","spat_window","estimate","pvalue","adj.r.squared")
  }else{
    estimates<- as.data.frame(matrix(nrow = (length(t_windows) * length(spat_windows)),ncol = 4))
    names(estimates) <- c("t_window","spat_window","estimate","pvalue")
  }
  matching <- as.data.frame(matrix(nrow = (length(t_windows) * length(spat_windows)),ncol = 10))
  names(matching) <- c("t_window","spat_window","control_pre","treatment_pre","L1_pre","commonSupport_pre","control_post","treatment_post","L1_post","commonSupport_post")
  SUTVA <- as.data.frame(matrix(nrow=0,ncol=8))
  estimates$spat_window <- spat_windows
  matching$spat_window <- spat_windows
  tmp <- c()
  for (i in 1:length(t_windows)){
    tmp <- c(tmp,rep(t_windows[i],length(spat_windows)))
  }
  estimates$t_window <- tmp
  matching$t_window <- tmp
  if (TCM){
    matchColumns <- c(matchColumns,"dependent_trend","SO_pre","MO_pre")
  } else {
    matchColumns <- c(matchColumns,"dependent_trend")
  }
  cem_cols <- c(matchColumns, "treatment", "dependent_post", "dependent_pre", estimationControls)
  for (space in seq(spat_windows[1],spat_windows[length(spat_windows)],spat_windows[2] - spat_windows[1])){
    for (time in seq(t_windows[1],t_windows[length(t_windows)],t_windows[2] - t_windows[1])){
      match_data <- subset(wakes,t_window == time & spat_window == space)
      matching$control_pre[matching$t_window == time & matching$spat_window == space] <- length(which(match_data$treatment == 0)) 
      matching$treatment_pre[matching$t_window == time & matching$spat_window == space] <- length(which(match_data$treatment == 1))
      matching$L1_pre[matching$t_window == time & matching$spat_window == space] <- round(imbalance(group = match_data$treatment, data = match_data[matchColumns])$L1[[1]],3)
      matching$commonSupport_pre[matching$t_window == time & matching$spat_window == space] <- round(imbalance(group = match_data$treatment, data = match_data[matchColumns])$L1[[3]],1)
      X <- match_data[,is.element(names(match_data),cem_cols)]
      if (match.default){
        mat <- do.call(cem,c(list(treatment = "treatment", data = X, drop = c("dependent_post", "dependent_pre", estimationControls)), ceminputs))
      }else{
        call <- call("cem.match")
        strata <- as.numeric(rep(NA,nrow(X)))
        n.strata <- as.integer(NA)
        vars <- matchColumns[!(matchColumns %in% c("treatment", "dependent_post", "dependent_pre", estimationControls))]
        drop <- c("treatment", "dependent_post", "dependent_pre", estimationControls)
        breaks <- list(match1 = NA, match2 = NA)
        treatment <- treatment
        n <- nrow(X)
        groups <- rep(0,nrow(X))
        groups[X$treatment=="1"] <- 1
        groups <- as.factor(groups)
        g.names <- c("0","1")
        n.groups <- 2
        G0 <- 1:nrow(X)
        G0 <- G0[X$treatment=="0"]
        G1 <- 1:nrow(X)
        G1 <- G1[X$treatment=="1"]
        group.idx <- list(G0 = G0, G1 = G1)
        group.len <- c(length(G0),length(G1))
        names(group.len) <- c("G0","G1")
        mstrata <- as.integer(NA,nrow(X))
        mstrataID <- NA
        matched <- rep(TRUE,nrow(X))
        baseline.group <- "1"
        tab <- as.table(rbind(c(length(G0),length(G1)),c(0,0),c(length(G0),length(G1))))
        colnames(tab) <- c("G0","G1")
        rownames(tab) <- c("All","Matched","Unmatched")
        k2k <- FALSE
        w <- rep(1,nrow(X))
        mat <- list(call = call, strata = strata, n.strata = n.strata, vars = vars, drop = drop, treatment = treatment, n = n, groups = groups, g.names = g.names, n.groups = n.groups, group.idx = group.idx, group.len = group.len, mstrata = mstrata, mstrataID = mstrataID, matched = matched, baseline.group = baseline.group, tab = tab, k2k = k2k, w = w)
        class(mat) <- "cem.match"
      }
      if (estimation == "lm"){
        if (weighted){
          md <- match_data[mat$matched,]
          weights <- mat$w[mat$matched]
          fit <- do.call(lm,c(list(formula = as.formula(formula), data = md, weights = weights), lminputs))
        }else{
          md <- match_data[mat$matched,]
          fit <- do.call(lm,c(list(formula = as.formula(formula), data = md), lminputs))
        }
      }
      if (estimation == "att"){
        if (weighted){
          md <- match_data[mat$matched,]
          fit <- do.call(att,c(list(obj = mat, formula = as.formula(formula), data = X), attinputs))
        }else{
          md <- match_data[mat$matched,]
          mat$w[mat$matched] <- rep(1,nrow(md))
          fit <- do.call(att,c(list(obj = mat, formula = as.formula(formula), data = X), attinputs))
        }
      }
      if (estimation == "nb"){
        md <- match_data[mat$matched,]
        fit <- do.call(glm.nb,c(list(formula = as.formula(formula), data = md, link = "identity"), nbinputs))
      }
      matching$control_post[matching$t_window == time & matching$spat_window == space] <- length(which(md$treatment == 0)) 
      matching$treatment_post[matching$t_window == time & matching$spat_window == space] <- length(which(md$treatment == 1))
      matching$L1_post[matching$t_window == time & matching$spat_window == space] <- round(imbalance(group = md$treatment, data = md[matchColumns])$L1[[1]],3)
      matching$commonSupport_post[matching$t_window == time & matching$spat_window == space] <- round(imbalance(group = md$treatment, data = md[matchColumns])$L1[[3]],1)
      if (estimation == "lm"){
        if (length(summary(fit)$coefficients[,1]) == nterms + 1){
          p_val <- summary(fit)$coefficients[nterms + 1,4]
          e_val <- summary(fit)$coefficients[nterms + 1,1]
          adj.r.squared <- summary(fit)$adj.r.squared
        } else {
          p_val <- 1
          e_val <- 0
          adj.r.squared <- 0
        }
      }
      if (estimation == "att"){
        if (length(fit$att.model[1,]) == nterms + 1){
          p_val <- fit$att.model[4,nterms + 1]
          e_val <- fit$att.model[1,nterms + 1]
        } else {
          p_val <- 1
          e_val <- 0
        }
      }
      if (estimation == "nb"){
        summary <- summary.glm(fit)
        if (length(summary$coefficients[,1]) == nterms + 1){
          p_val <- summary$coefficients[nterms+1,4]
          e_val <- summary$coefficients[nterms+1,1]
        } else {
          p_val <- 1
          e_val <- 0
        }
      }  
      if (is.na(p_val) || is.na(e_val) || is.nan(p_val) || is.nan(e_val)){
        p_val <- 1
        e_val <- 0
      }
      estimates$estimate[estimates$t_window == time & estimates$spat_window == space] <- e_val
      estimates$pvalue[estimates$t_window == time & estimates$spat_window == space] <- p_val
      if (estimation == "lm"){
        estimates$adj.r.squared[estimates$t_window == time & estimates$spat_window == space] <- adj.r.squared
      }
      sub <- subset(wakes,(wakes$t_window==time & wakes$spat_window==space))
      doubles_pre_sub <- subset(wakes, (wakes$t_window==time & wakes$spat_window==space & wakes$SO_pre>0))
      doubles_pre <- round(nrow(doubles_pre_sub)/nrow(sub),3)
      doubles_post_sub <- subset(wakes, (wakes$t_window==time & wakes$spat_window==space & wakes$SO_post>0))
      doubles_post <- round(nrow(doubles_post_sub)/nrow(sub),3)
      doubles_sub <- subset(wakes, (wakes$t_window==time & wakes$spat_window==space & (wakes$SO_pre>0 | wakes$SO_post>0)))
      doubles <- round(nrow(doubles_sub)/nrow(sub),3)
      spills_pre_sub <- subset(wakes, (wakes$t_window==time & wakes$spat_window==space & wakes$MO_pre>0))
      spills_pre <- round(nrow(spills_pre_sub)/nrow(sub),3)
      spills_post_sub <- subset(wakes, (wakes$t_window==time & wakes$spat_window==space & wakes$MO_post>0))                                            
      spills_post <- round(nrow(spills_post_sub)/nrow(sub),3)
      spills_sub <- subset(wakes, (wakes$t_window==time & wakes$spat_window==space & (wakes$MO_pre>0 | wakes$MO_post>0)))                                              
      spills <- round(nrow(spills_sub)/nrow(sub),3)
      SUTVA <- rbind(SUTVA,c(time,space,doubles_pre,doubles_post,doubles,spills_pre,spills_post,spills))
    }
  }
  names(SUTVA) <- c("t_window","spat_window","SO_pre","SO_post","SO","MO_pre","MO_post","MO")
  SUTVA <- SUTVA[order(SUTVA$t_window),]
  row.names(SUTVA) <- NULL
  matchedwake <- list (estimates = estimates, matching = matching, SUTVA = SUTVA, wakes = wakes)
  return(matchedwake)  
}

slidingWake <- function(data, t_unit, t_window, spat_window, treatment, control, dependent, matchColumns, estimationControls, memory){
  data <- data[order(data$timestamp),] 
  row.names(data)<-NULL
  datatype<-as.character(lapply(1:length(data),function(x) class(data[,x])))
  treatmentinput<-treatment
  treatmentinput[1]<-as.character(grep(treatment[1],names(data)))
  controlinput<-control
  controlinput[1]<-as.character(grep(control[1],names(data)))
  depvarinput<-dependent
  depvarinput[1]<-as.character(grep(dependent[1],names(data)))
  wakeparams<-array(0,4)
  wakeparams[1]<-grep('lat',names(data))
  wakeparams[2]<-grep('long',names(data))
  wakeparams[3]<-grep('timestamp',names(data))
  wakeparams<-as.character(wakeparams)
  wakeparams[4]<-t_unit
  datainput<-as.matrix(data)
  mode(datainput)<-'character'
  if (t_unit=='days'){
    data$timestamp <- strptime(data$timestamp, "%Y-%m-%d %H:%M:%S")
    data$timestamp <- format(data$timestamp,"%Y-%m-%d 00:00:00")
  }
  if (t_unit=='hours'){
    data$timestamp <- strptime(data$timestamp, "%Y-%m-%d %H:%M:%S")
    data$timestamp <- format(data$timestamp,"%Y-%m-%d %H:00:00")
  }
  if (t_unit=='mins'){
    data$timestamp <- strptime(data$timestamp, "%Y-%m-%d %H:%M:%S")
    data$timestamp <- format(data$timestamp,"%Y-%m-%d %H:%M:00")
  }
  wake_events<-subset(data, data[,c(treatment[1])]==treatment[2])
  wake_events<-rbind(wake_events,subset(data, data[,c(control[1])]==control[2]))
  wake_events[wake_events[,c(treatment[1])]==treatment[2],c(treatmentinput[1])]<-'1'
  wake_events[wake_events[,c(control[1])]==control[2],c(treatmentinput[1])]<-'0'
  wakeindex<-cbind(row.names(wake_events),wake_events[,c(treatmentinput[1])])
  wakeindex<- wakeindex[order(as.numeric(wakeindex[,1])), ]
  timevarinput <- as.character(seq(t_window[1],t_window[2],length=(t_window[2]-t_window[1])/t_window[3]+1))
  spatvarinput <- as.character(seq(spat_window[1],spat_window[2],length=(spat_window[2]-spat_window[1])/spat_window[3]+1))
  if (length(estimationControls) > 0){
    matchColums <- c(matchColumns,estimationControls)
  }
  matchCol <- as.character(lapply(matchColumns,function(x) grep(x,names(data))))
  directory <- getwd()
  cat(paste('Initializing JVM (',memory,'GB memory)...',sep=""))
  .jinit(parameters=paste("-Xmx",memory,"g",sep=""),force.init = TRUE)
  javatimevarinput <- .jcast(.jarray(timevarinput,contents.class="[Ljava/lang/String;",dispatch=TRUE),"[Ljava/lang/String;")      
  javaspatvarinput <- .jcast(.jarray(spatvarinput,contents.class="[Ljava/lang/String;",dispatch=TRUE),"[Ljava/lang/String;")      
  javamatchCol <- .jcast(.jarray(matchCol,contents.class="[Ljava/lang/String;",dispatch=TRUE),"[Ljava/lang/String;")      
  javadata<-.jcast(.jarray(datainput,contents.class="[[Ljava/lang/String;",dispatch=TRUE),"[[Ljava/lang/String;")
  javawakeindex<-.jcast(.jarray(wakeindex,contents.class="[[Ljava/lang/String;",dispatch=TRUE),"[[Ljava/lang/String;")  
  javadirectory <- .jcast(.jarray(directory,contents.class="[Ljava/lang/String;",dispatch=TRUE),"[Ljava/lang/String;")      
  obj <- .jnew("WakeCounter")
  if(memory<10){
    cat('                                          [OK]\n')
  }else{
    cat('                                         [OK]\n')
  }
  cat('Iterating through sliding windows and generating wakes...')
  out <- .jcall(obj,"[[Ljava/lang/String;","wakeCounting", javatimevarinput, javaspatvarinput, javadata, javawakeindex, wakeparams, treatmentinput, controlinput, depvarinput, javamatchCol, javadirectory)
  wakes <- lapply(1:length(out),function(x) .jevalArray(out[[x]],simplify=TRUE))
  wakes<-do.call(rbind, wakes)
  wakes<-as.data.frame(wakes)
  names(wakes) <- c("eventID","t_window", "spat_window", "treatment","dependent_pre", "dependent_trend", "SO_pre", "MO_pre", "dependent_post", "SO_post", "MO_post", matchColumns)
  wakesoutput<-subset(wakes,wakes$t_window!="NaN")
  row.names(wakesoutput)<-NULL
  cat('                 [OK]\n')
  if (length(wakesoutput[,1])!=length(wakes[,1])){
    cat('WARNING: Some wakes were incomplete and had to be disregarded; please see matchedwake_log.txt for details.\n')  
  }
  dataformat<-c("numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric")
  matchCol<-as.numeric(matchCol)
  dataformat<-c(dataformat,datatype[matchCol])
  for (row in 1:length(dataformat)){
    if (dataformat[row]=="numeric" || dataformat[row]=="integer"){
      wakesoutput[,row] <- as.numeric(as.character(wakesoutput[,row]))
    }
    if (dataformat[row]=="factor"){
      wakesoutput[,row] <- as.factor(as.character(wakesoutput[,row]))
    }
    if (dataformat[row]=="character"){
      wakesoutput[,row] <- as.character(wakesoutput[,row])
    }
  }
  return(wakesoutput)
}
