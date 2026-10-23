# GMM provenance and supplementary diagnostics. Source without side effects.
gd_rng <- function() if (exists('.Random.seed', .GlobalEnv, inherits=FALSE))
  get('.Random.seed', .GlobalEnv) else NULL

gd_seed <- function(seed, expr) {
  old <- gd_rng()
  on.exit(if (is.null(old)) {
    if (exists('.Random.seed', .GlobalEnv, inherits=FALSE)) rm('.Random.seed', envir=.GlobalEnv)
  } else assign('.Random.seed', old, envir=.GlobalEnv))
  set.seed(seed)
  force(expr)
}

gd_draw <- function(pool, n, seed) gd_seed(seed, pool[sample.int(length(pool), min(n,length(pool)))])
gd_csv <- function(x, path) utils::write.csv(x, path, row.names=FALSE, na='')

gd_capture <- function(expr) {
  warnings <- character(); error <- NA_character_
  fit <- tryCatch(withCallingHandlers(force(expr), warning=function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
  }), error=function(e) { error <<- conditionMessage(e); NULL })
  if(is.null(fit) && is.na(error)) error <- 'Mclust returned NULL (no fitted model; possible singularity)'
  list(fit=fit, warnings=warnings, error=error)
}

gd_compact <- function(fit) {
  if (is.null(fit)) return(NULL)
  fit[c('data','z','classification','uncertainty')] <- NULL
  fit
}

gd_model_ok <- function(model) {
  if (is.null(model) || !is.finite(model$loglik) || any(!is.finite(model$parameters$mean))) return(FALSE)
  s <- model$parameters$variance$sigma
  if (is.null(s) || length(dim(s)) != 3 || any(!is.finite(s))) return(FALSE)
  all(vapply(seq_len(dim(s)[3]), function(k)
    isTRUE(all.equal(s[,,k],t(s[,,k]),tolerance=1e-10)) &&
      min(eigen(s[,,k],symmetric=TRUE,only.values=TRUE)$values)>0, logical(1))) &&
    all(model$parameters$pro>0) && abs(sum(model$parameters$pro)-1)<1e-8
}

gd_validate_z <- function(z) {
  if (!is.matrix(z) || ncol(z)<2 || any(!is.finite(z)) || any(z<0) || any(z>1) ||
      any(abs(rowSums(z)-1)>1e-8)) stop('Invalid posterior probability matrix')
  invisible(TRUE)
}

gd_z <- function(model, X) {
  if (!gd_model_ok(model)) stop('Invalid fitted Gaussian parameters')
  z <- predict(model, X)$z
  gd_validate_z(z)
  z
}

gd_perms <- function(x) {
  if(length(x)==1) return(matrix(x, nrow=1))
  do.call(rbind,lapply(seq_along(x),function(i) cbind(x[i],gd_perms(x[-i]))))
}

gd_compare <- function(reference, candidate) {
  gd_validate_z(reference); gd_validate_z(candidate)
  stopifnot(nrow(reference)==nrow(candidate))
  r <- max.col(reference,ties.method='first'); c <- max.col(candidate,ties.method='first')
  ari <- mclust::adjustedRandIndex(r,c)
  if (!is.finite(ari)) {
    contingency <- table(r,c)>0
    if(all(rowSums(contingency)==1) && all(colSums(contingency)==1)) ari <- 1
  }
  if(ncol(reference)!=ncol(candidate)) return(list(ari=ari,agreement=NA_real_,posterior_mae=NA_real_,order=NULL))
  perms <- gd_perms(seq_len(ncol(reference)))
  scores <- apply(perms,1,function(p) sum(r==match(c,p)))
  ord <- perms[which.max(scores),]
  aligned <- candidate[,ord,drop=FALSE]
  list(ari=ari, agreement=mean(r==max.col(aligned,ties.method='first')),
       posterior_mae=mean(abs(reference-aligned)),order=ord)
}

gd_archive <- function(directory, model, candidates, rows, normalisation,
                       rng, metadata, keys=NULL) {
  # Do not draw random numbers here: observation must not alter later pipeline RNG.
  if(dir.exists(directory)) stop('Archive directory already exists: ',directory)
  dir.create(directory,recursive=TRUE)
  payload <- list(schema_version=1L, model=model, candidates=candidates,
                  rows=rows, keys=keys, normalisation=normalisation, rng=rng,
                  metadata=metadata, session=utils::sessionInfo(), rng_kind=RNGkind())
  saveRDS(payload,file.path(directory,'gmm_archive.rds'),compress='gzip')
  writeLines(capture.output(utils::sessionInfo()),file.path(directory,'sessionInfo.txt'))
  gd_csv(data.frame(file='gmm_archive.rds',md5=unname(tools::md5sum(file.path(directory,'gmm_archive.rds')))),
         file.path(directory,'archive_checksum.csv'))
  invisible(directory)
}

gd_fit <- function(X, G, seed) gd_seed(seed, {
  before <- gd_rng()
  result <- gd_capture(mclust::Mclust(X,G=G,modelNames='VVV',verbose=FALSE))
  result$seed <- seed; result$rng_before <- before; result$rng_after <- gd_rng()
  result
})

gd_candidate_row <- function(x,G,n,seed) data.frame(G=G,n=n,seed=seed,
  success=gd_model_ok(x$fit), BIC=if(is.null(x$fit)) NA_real_ else as.numeric(x$fit$bic),
  loglik=if(is.null(x$fit)) NA_real_ else x$fit$loglik,
  df=if(is.null(x$fit)) NA_integer_ else x$fit$df,
  warnings=paste(x$warnings,collapse=' | '), error=x$error)

gd_selection <- function(X,rows,seeds,directory) {
  candidates <- list(); tab <- data.frame()
  for(g in 2:5) {
    result <- gd_fit(X[rows,,drop=FALSE],g,seeds[g])
    candidates[[as.character(g)]] <- list(model=gd_compact(result$fit),warnings=result$warnings,
      error=result$error,seed=seeds[g],rng_before=result$rng_before,rng_after=result$rng_after)
    tab <- rbind(tab,gd_candidate_row(result,g,length(rows),seeds[g]))
    gd_csv(tab,file.path(directory,'candidate_models.csv'))
    saveRDS(candidates,file.path(directory,'candidate_models.rds'))
    rm(result); gc(FALSE)
  }
  eligible <- tab$success & is.finite(tab$BIC)
  if(!all(eligible[tab$G %in% 3:5])) stop('At least one canonical candidate failed; selection incomplete')
  tab$delta_BIC_from_expanded_best <- max(tab$BIC[eligible])-tab$BIC
  canonical <- tab$G %in% 3:5
  tab$delta_BIC_from_canonical_best <- max(tab$BIC[eligible&canonical])-tab$BIC
  gd_csv(tab,file.path(directory,'candidate_models.csv'))
  list(canonical_G=tab$G[which.max(replace(tab$BIC,!canonical,-Inf))],
       expanded_G=tab$G[which.max(replace(tab$BIC,!eligible,-Inf))],candidates=candidates)
}

gd_uncertainty <- function(z,keys,directory,prefix='saved') {
  gd_validate_z(z)
  largest <- do.call(pmax,as.data.frame(z))
  second <- if(ncol(z)==3) rowSums(z)-largest-do.call(pmin,as.data.frame(z)) else
    apply(z,1,function(x) sort(x,decreasing=TRUE)[2])
  entropy <- -rowSums(z*log(pmax(z,.Machine$double.xmin)))/log(ncol(z))
  describe <- function(i,dimension,group) data.frame(dimension=dimension,group=group,n=length(i),
    mean_uncertainty=mean(1-largest[i]),mean_normalised_entropy=mean(entropy[i]),
    median_margin=median(largest[i]-second[i]),
    share_max_below_0.5=mean(largest[i]<.5),share_max_below_0.8=mean(largest[i]<.8),
    share_max_below_0.9=mean(largest[i]<.9))
  result <- describe(seq_len(nrow(z)),'all','all')
  for(d in c('vessel','year')) for(g in unique(keys[[d]])) {
    result <- rbind(result,describe(which(keys[[d]]==g),d,as.character(g)))
  }
  gd_csv(result,file.path(directory,paste0(prefix,'_uncertainty.csv')))
}

gd_density_exports <- function(model,X,directory) {
  stopifnot(gd_model_ok(model))
  m <- model$parameters$mean; s <- model$parameters$variance$sigma; w <- model$parameters$pro
  ord <- order(m[1,]); g <- length(w)
  means <- covariances <- numerical <- tails <- curves <- bins <- ellipses <- data.frame()
  limits <- list(c(.5,20),c(0,1),c(0,1))
  for(k in seq_len(g)) {
    j <- ord[k]; eig <- eigen(s[,,j],symmetric=TRUE,only.values=TRUE)$values
    numerical <- rbind(numerical,data.frame(component=k,original_component=j,weight=w[j],
      eigen_min=min(eig),eigen_max=max(eig),condition_number=max(eig)/min(eig)))
    for(v in 1:3) means <- rbind(means,data.frame(component=k,variable=colnames(X)[v],mean=m[v,j],sd=sqrt(s[v,v,j])))
    for(v in 1:3) for(u in 1:3) covariances <- rbind(covariances,
      data.frame(component=k,variable1=colnames(X)[v],variable2=colnames(X)[u],covariance=s[v,u,j]))
    for(pair in list(c(1,2),c(1,3),c(2,3))) {
      eg <- eigen(s[pair,pair,j],symmetric=TRUE)
      angle <- seq(0,2*pi,length.out=181)
      circle <- rbind(cos(angle),sin(angle))
      for(level in c(.5,.9)) {
        xy <- eg$vectors %*% diag(sqrt(eg$values*qchisq(level,df=2)),2) %*% circle
        xy <- sweep(xy,1,m[pair,j],'+')
        ellipses <- rbind(ellipses,data.frame(component=k,variable1=colnames(X)[pair[1]],
          variable2=colnames(X)[pair[2]],probability_content=level,vertex=seq_along(angle),
          x=xy[1,],y=xy[2,]))
      }
    }
  }
  for(v in 1:3) {
    lo <- limits[[v]][1]; hi <- limits[[v]][2]; sdv <- sqrt(s[v,v,])
    xx <- seq(min(lo,m[v,]-4*sdv),max(hi,m[v,]+4*sdv),length.out=1001)
    dd <- sapply(seq_len(g),function(k) w[k]*dnorm(xx,m[v,k],sdv[k]))
    for(k in seq_len(g)) curves <- rbind(curves,data.frame(variable=colnames(X)[v],x=xx,
      component=k,density_contribution=dd[,ord[k]],mixture_density=rowSums(dd)))
    for(k in seq_len(g)) {
      j <- ord[k]
      tails <- rbind(tails,data.frame(variable=colnames(X)[v],component=k,weight=w[j],
        lower_bound=lo,upper_bound=hi,prob_below=pnorm(lo,m[v,j],sdv[j]),
        prob_above=pnorm(hi,m[v,j],sdv[j],lower.tail=FALSE),
        observed_at_lower=mean(X[,v]==lo),observed_at_upper=mean(X[,v]==hi)))
    }
    breaks <- seq(lo,hi,length.out=101)
    counts <- hist(X[,v],breaks=breaks,plot=FALSE,include.lowest=TRUE)$counts
    cdf <- function(q) as.numeric(sapply(seq_len(g),function(k) w[k]*pnorm(q,m[v,k],sdv[k])) %*% rep(1,g))
    empirical <- ecdf(X[,v])
    bins <- rbind(bins,data.frame(variable=colnames(X)[v],left=head(breaks,-1),right=tail(breaks,-1),
      observed_probability=counts/nrow(X),model_probability=diff(cdf(breaks)),
      empirical_cdf=empirical(tail(breaks,-1)),model_cdf=cdf(tail(breaks,-1))))
  }
  for(n in c('means','covariances','numerical','tails','curves','bins','ellipses'))
    gd_csv(get(n),file.path(directory,paste0('model_',n,'.csv')))
}

gd_partition <- function(keys,dimension,group,train_n,eval_n,seed) {
  test <- gd_draw(which(keys[[dimension]]==group),eval_n,seed)
  train <- gd_draw(which(keys[[dimension]]!=group),train_n,seed+1L)
  control <- gd_draw(setdiff(seq_len(nrow(keys)),test),length(train),seed+2L)
  if(length(intersect(train,test)) || any(keys[[dimension]][train]==group) ||
     length(intersect(control,test))) stop('Holdout isolation failed')
  list(train=train,test=test,control=control)
}

gd_stability <- function(X,keys,model,directory,seed=20260924L,reps=10L,train_n=100000L,eval_n=100000L) {
  log <- data.frame(); manifest <- list()
  append_result <- function(name,kind,group,fit,rows,eval,reference,against) {
    path <- paste0('stability_',name,'.rds')
    saveRDS(list(model=gd_compact(fit$fit),warnings=fit$warnings,error=fit$error,
                 seed=fit$seed,rng_before=fit$rng_before,rng_after=fit$rng_after,
                 train_rows=rows,evaluation_rows=eval),file.path(directory,path))
    valid <- gd_model_ok(fit$fit)
    cmp <- if(valid) gd_compare(reference,gd_z(fit$fit,X[eval,,drop=FALSE])) else
      list(ari=NA_real_,agreement=NA_real_,posterior_mae=NA_real_)
    log <<- rbind(log,data.frame(name=name,kind=kind,group=group,against=against,
       train_n=length(rows),evaluation_n=length(eval),success=valid,ARI=cmp$ari,
       agreement=cmp$agreement,posterior_mae=cmp$posterior_mae,
       component_order=if(is.null(cmp$order)) NA_character_ else paste(cmp$order,collapse=','),
       seed=fit$seed,warnings=paste(fit$warnings,collapse=' | '),error=fit$error,archive=path))
    gd_csv(log,file.path(directory,'stability.csv'))
  }
  evaluation <- gd_draw(seq_len(nrow(X)),min(eval_n,floor(nrow(X)/5)),seed)
  pool <- setdiff(seq_len(nrow(X)),evaluation)
  fixed_rows <- gd_draw(pool,train_n,seed+1L)
  anchor <- gd_fit(X[fixed_rows,,drop=FALSE],model$G,seed+2L)
  if(!gd_model_ok(anchor$fit)) {
    append_result('anchor','anchor','all',anchor,fixed_rows,evaluation,
                  gd_z(model,X[evaluation,,drop=FALSE]),'canonical')
    stop('Stability anchor failed')
  }
  anchor_z <- gd_z(anchor$fit,X[evaluation,,drop=FALSE])
  append_result('anchor','anchor','all',anchor,fixed_rows,evaluation,gd_z(model,X[evaluation,,drop=FALSE]),'canonical')
  for(r in seq_len(reps)) {
    init <- gd_fit(X[fixed_rows,,drop=FALSE],model$G,seed+10L+r)
    append_result(paste0('init_',r),'initialisation','all',init,fixed_rows,evaluation,anchor_z,'same_size_anchor')
    new_rows <- gd_draw(pool,train_n,seed+100L+r)
    samp <- gd_fit(X[new_rows,,drop=FALSE],model$G,seed+2L)
    append_result(paste0('sample_',r),'resampling','all',samp,new_rows,evaluation,anchor_z,'same_size_anchor')
  }
  count <- 0L
  for(d in c('vessel','year')) for(g in unique(keys[[d]])) {
    count <- count+1L; ss <- seed+1000L+10L*count
    part <- gd_partition(keys,d,g,train_n,min(eval_n,20000L),ss)
    manifest[[paste(d,g,sep='_')]] <- c(part,list(seed=ss))
    ref <- gd_z(model,X[part$test,,drop=FALSE])
    f <- gd_fit(X[part$train,,drop=FALSE],model$G,ss+3L)
    ctrl <- gd_fit(X[part$control,,drop=FALSE],model$G,ss+3L)
    append_result(paste0(d,'_',count),'group_exclusion',as.character(g),f,part$train,part$test,ref,'canonical')
    append_result(paste0(d,'_control_',count),'all_group_control',as.character(g),ctrl,part$control,part$test,ref,'canonical')
    if(gd_model_ok(ctrl$fit)) append_result(paste0(d,'_paired_',count),'group_exclusion',as.character(g),
        f,part$train,part$test,gd_z(ctrl$fit,X[part$test,,drop=FALSE]),'same_size_all_group_control')
    saveRDS(manifest,file.path(directory,'group_partitions.rds'))
  }
  invisible(log)
}
