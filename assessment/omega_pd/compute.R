# Omega x excavation depth. Execute ONLY after approval of the rule (RULE.md).
# Equations: authoritative merge/step6, full md5 frozen in inputs.tsv.
# Numeric indexing, logical counts, and exit codes are programming constants;
# all scientific constants/thresholds below carry their provenance at assignment.
suppressPackageStartupMessages({library(arrow); library(data.table); library(yaml)})
T <- path.expand(file.path(Sys.getenv("DFP_REPO"), "assessment/omega_pd"))
OUT <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_incert/omega_pd")
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Slurm allocation required")
if (!dir.exists(file.path(OUT, "run_claim"))) stop("Run claim required")
ZERO <- 0 # Algebraic zero, as in the merge and step6.
ONE <- 1 # Unit fractions/caps, as in the merge and step6. Also one year.
HALF_LAYER_M <- 0.05 # Thickness of each upper slice, metres (as in the merge).
MIXED_M <- 0.10 # Mixed layer bottom, metres (as in the merge).
CM_PER_M <- 100 # SI unit conversion, as in the merge.
M2_PER_KM2 <- 1e6 # Exact SI area conversion, for table cell_km2 and V_balaye in m3.
PD_REF <- 1 # Central depth, m; checked cellwise at runtime.
PD <- c(0.3, 0.5, 1) # Minimum of the grid; 0.3 from Porz et al. 2026 (Methods); 0.5 design midpoint; 1 central.
OMEGA <- c(0.1, 0.42, 1, 5) # Omega scenarios of the YAML in service. Trials, not measured site values.
OMEGA_REF <- 0.42 # Central value of the YAML; Middelburg attribution only via the YAML, not rechecked here.
OMEGA_HIGH <- 5 # YAML scenario, trial used for the primary contrast; not an empirical port value.
REL <- 1e-9 # Design choice, no source: per-cell relative numerical tolerance.
ABS <- 1e-12 # Design choice, no source: absolute floor for near-zero cell values.
TOTAL_REL <- 1e-6 # Lock V0 of the uncertainty tests (read.R); reused for totals and unrounded Omega ratios.
ROUND_TOL <- 5.1e-5 # Design choice, no source: half last displayed decimal + margin.
MAT <- 0.05 # 5 % practical convention of the uncertainty tests, extended to interaction by design.
PRINT_DIGITS <- 15 # Design choice: reproducible reporting precision, not a test threshold.
options(digits = PRINT_DIGITS)
LABELS <- c("smooth", "flag")
ARMS <- c("frac_fixe", "frac_recalculee")
CASE_OMEGA <- c(ref=OMEGA_REF, "om0.1"=OMEGA[1], "om1.0"=OMEGA[3], "om5.0"=OMEGA_HIGH) # Scenario order from the YAML; indices only.
EXPECTED <- data.table(
  label=rep(LABELS, each=length(OMEGA)), omega=rep(OMEGA, times=length(LABELS)),
  ratio=c(0.8958, 1, 1.1872, 1.9340, 0.8958, 1, 1.1870, 1.9335)) # Omega ratios of the uncertainty tests (read.R), rounded to four decimals.
TOTAL_EXPECTED <- c(smooth=4038.1542, flag=5654.3864) # Reference totals (t C); four decimals, tested with ROUND_TOL.

write_new <- function(lines, file) {
  target <- file.path(OUT, file)
  if (file.exists(target)) stop("Output already exists: ", file)
  writeLines(lines, target, useBytes=TRUE)
}
save_table <- function(x, file) {
  target <- file.path(OUT, file)
  if (file.exists(target)) stop("Output already exists: ", file)
  fwrite(x, target, sep="\t", na="NA")
}
fail <- function(message) stop(message, call.=FALSE)
need <- function(x, columns, name) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) fail(paste("Missing columns in", name, paste(absent, collapse=", ")))
  if (!nrow(x) || anyNA(x$grid_id) || anyDuplicated(x$grid_id)) fail(paste("Empty/invalid/duplicate grid_id:", name))
}
finite <- function(x, name, lower=ZERO, upper=Inf, allow_na=FALSE) {
  use <- if (allow_na) !is.na(x) else rep(TRUE, length(x))
  if ((!allow_na && anyNA(x)) || any(!is.finite(x[use])) || any(x[use] < lower | x[use] > upper)) fail(paste("Invalid values:", name))
}
close <- function(a, b, name, rel=REL, abs=ABS) {
  if (length(a) != length(b) || any(is.na(a) != is.na(b))) fail(paste("Length/NA mismatch:", name))
  ii <- !is.na(a)
  if (any(!is.finite(a[ii])) || any(!is.finite(b[ii]))) fail(paste("Nonfinite:", name))
  delta <- base::abs(a[ii] - b[ii])
  if (any(delta > abs + rel * base::abs(b[ii]))) fail(paste("Numerical mismatch:", name))
  if (length(delta)) max(delta) else ZERO
}
align <- function(x, ids, name) {
  if (nrow(x) != length(ids) || !setequal(x$grid_id, ids)) fail(paste("Different grid_id universe:", name))
  x[match(ids, grid_id)]
}
locks <- list()
record <- function(label, check, error=ZERO) {
  locks[[length(locks)+1L]] <<- data.table(label=label, check=check, max_error=error, status="PASS")
}

main <- function() {
  inputs <- fread(file.path(T, "inputs.tsv"), colClasses="character")
  if (anyDuplicated(inputs$id)) fail("Duplicated input identifier")
  digest <- unname(tools::md5sum(inputs$path))
  if (anyNA(digest) || any(digest != inputs$md5)) fail("Input fingerprint changed")
  file_of <- function(id) {
    ans <- inputs$path[match(id, inputs$id)]
    if (is.na(ans)) fail(paste("Unknown input", id))
    ans
  }
  read_pq <- function(id) as.data.table(arrow::read_parquet(file_of(id)))
  pars <- yaml::read_yaml(file_of("yaml"))$scenarios$default
  if (!identical(pars$depth_mode, "mechanistic")) fail("Wrong depth mode")
  close(pars$omega_cm_yr, OMEGA_REF, "reference omega")
  # Extract only the two top-level function definitions; never source the production script.
  expressions <- parse(file=file_of("merge"), keep.source=TRUE)
  env <- new.env(parent=baseenv())
  for (nm in c("weight_depth", "weight_eff")) {
    found <- Filter(function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
      identical(e[[2L]], as.name(nm)), as.list(expressions))
    if (length(found) != ONE) fail(paste("Function missing/ambiguous:", nm))
    eval(found[[1L]], envir=env) # Sole verified definition, no execution of pipeline main.
  }
  env$Db_cm2 <- pars$Db_cm2_yr # Held fixed across cases (YAML).
  env$fast_frac <- pars$fast_fraction # Held fixed (YAML).
  env$omega_cm <- pars$omega_cm_yr
  finite(c(pars$Db_cm2_yr, pars$slow_k, pars$preservation_factor), "YAML constants")
  if (any(c(pars$Db_cm2_yr, pars$slow_k) <= ZERO)) fail("Positive rates required")
  finite(pars$fast_fraction, "pool fraction", upper=ONE)
  record("all", "source_fingerprints_and_functions")

  calc <- function(d, omega, pd, arm) {
    env$omega_cm <- omega
    # Production passes one depth per cell. In weight_depth, ifelse(H3>0,...)
    # takes its output length from H3: a scalar would truncate the k-dependent deep term.
    pd <- rep(pd, nrow(d))
    # Same operations and order as the merge; constants sourced above.
    w1 <- pmin(HALF_LAYER_M, pd) / pd
    w2 <- pmin(HALF_LAYER_M, pmax(pd - HALF_LAYER_M, ZERO)) / pd
    w3 <- pmax(pd - MIXED_M, ZERO) / pd
    H3 <- pmax(pd - MIXED_M, ZERO) * CM_PER_M
    W <- env$weight_eff(d$k_used, pars$slow_k, w1, w2, w3, H3, t=ONE)
    # Full reconstruction, not multiplication of a previously capped f_i.
    svr <- d$SAR * pd
    # In the merge, sum_dw=0 gives p_d=NA and SVR=NA. Preserve that NA mask;
    # the input gate below permits it only in zero-SAR cells, not in active cells.
    svr[is.na(d$SVR)] <- NA_real_
    frac <- d$frac_sur_place
    if (arm == "frac_recalculee") {
      matched <- !is.na(d$V_retire_test)
      swept <- svr[matched] * d$cell_km2_test[matched] * M2_PER_KM2
      if (any(swept <= ZERO)) fail("Nonpositive swept volume on matched hopper cells")
      frac[matched] <- pmax(swept - d$V_retire_test[matched], ZERO) / swept
      # Unmatched cells retain the production convention frac=1 (as in the merge).
    }
    pl_corr <- d$p_l * W * d$fresh_fact
    before <- svr * pl_corr * pars$preservation_factor *
      (pars$fast_fraction * (ONE-exp(-d$k_used*ONE)) +
       (ONE-pars$fast_fraction) * (ONE-exp(-pars$slow_k*ONE)))
    raw <- before * frac # Order BEFORE cap, as in the merge.
    fi_merge <- pmin(pmax(raw, ZERO), ONE)
    fi <- fi_merge
    fi[is.na(fi)] <- ZERO # As in step6. Counted explicitly in diagnostics.
    cpre <- d$C0i * fi * d$di # As in step6, di guarded at one.
    cri <- pmin(cpre, d$C0i)
    data.table(grid_id=d$grid_id, W_eff=W, SVR=svr, frac=frac, before_frac=before,
               raw_f=raw, fi_merge=fi_merge, fi=fi, C0i=d$C0i, cpre=cpre, C_ri=cri)
  }

  data <- list()
  for (lab in LABELS) {
    a <- read_pq(paste0("central_fi_", lab))
    # p_d_safe is not saved (the merge builds it, then drops it); rebuilt below with the merge's rule.
    need(a, c("grid_id","SAR","SVR","p_d","p_l","k_used","fresh_fact",
              "W_eff","f_i_full","f_i_full_avant_pr","frac_sur_place"), paste(lab,"fi"))
    b <- read_pq(paste0("central_cri_", lab))
    need(b, c("grid_id","C0i","C0i_mean_carbon_stock","di","f_i_full","C_ri"), paste(lab,"cri"))
    b <- align(b, a$grid_id, paste(lab,"central cri"))
    finite(a$SAR,"SAR"); finite(a$SVR,"SVR",allow_na=TRUE)
    finite(a$p_l,"p_l",allow_na=TRUE); finite(a$fresh_fact,"fresh_fact")
    finite(a$k_used,"k_used"); if(any(a$k_used <= ZERO)) fail("Nonpositive k")
    finite(b$C0i,"C0i"); finite(b$C_ri,"central C_ri")
    close(b$C0i,b$C0i_mean_carbon_stock,"mean stock alias")
    close(b$di, rep(ONE,nrow(b)), "di=1")
    active <- a$SAR > ZERO
    close(a$p_d[active],rep(PD_REF,sum(active)),"all active central depths are 1m")
    svr_ref <- a$SAR*PD_REF
    svr_ref[a$SAR==ZERO & is.na(a$p_d)] <- NA_real_ # Same zero-volume convention as the merge.
    close(a$SVR,svr_ref,"central SVR=SAR*1m with inactive NA convention")
    p_d_safe <- fifelse(is.na(a$p_d) | a$p_d <= ZERO, ONE, a$p_d) # As in the merge: NA/0 -> 1 m.
    close(p_d_safe,rep(PD_REF,nrow(a)),"central safe depth")
    a[, `:=`(C0i=b$C0i, di=b$di)]
    fr <- fread(file_of(paste0("frac_",lab)))
    # The complete schema of the table in service; no header is printed.
    need(fr,c("grid_id","key","SAR","SVR","cell_km2","V_balaye","V_retire","V_sur_place","frac"),paste(lab,"frac"))
    for (col in c("SAR","SVR","cell_km2","V_balaye","V_retire","V_sur_place")) finite(fr[[col]],col)
    finite(fr$frac,"frac",upper=ONE)
    if(any(fr$V_balaye <= ZERO) || any(fr$cell_km2 <= ZERO)) fail("Nonpositive source volume/area")
    close(fr$V_balaye,fr$SVR*fr$cell_km2*M2_PER_KM2,"table swept-volume units")
    close(fr$V_sur_place,pmax(fr$V_balaye-fr$V_retire,ZERO),"remaining volume")
    close(fr$frac,pmax(fr$V_balaye-fr$V_retire,ZERO)/fr$V_balaye,"source frac reconstruction")
    ix <- match(a$grid_id, fr$grid_id)
    matched <- !is.na(ix)
    close(a$SAR[matched],fr$SAR[ix[matched]],"table/central SAR")
    close(a$SVR[matched],fr$SVR[ix[matched]],"table/central SVR")
    expected_frac <- rep(ONE,nrow(a)) # Production default for absent grid_id, as in the merge.
    expected_frac[matched] <- fr$frac[ix[matched]]
    close(a$frac_sur_place,expected_frac,"central fraction with explicit absent-row default")
    a[, `:=`(V_retire_test=fr$V_retire[ix], cell_km2_test=fr$cell_km2[ix])]
    data[[lab]] <- a
    record(lab,"schemas_unique_ids_depth_stock_frac")
    for (arm in ARMS) {
      v <- calc(a,OMEGA_REF,PD_REF,arm)
      record(lab,paste(arm,"central_W"),close(v$W_eff,a$W_eff,"central W"))
      record(lab,paste(arm,"central_before_frac"),close(v$before_frac,a$f_i_full_avant_pr,"uncapped pre-fraction f"))
      record(lab,paste(arm,"central_fi"),close(v$fi_merge,a$f_i_full,"central fi"))
      close(v$fi,b$f_i_full,"step6 central fi")
      record(lab,paste(arm,"central_cri"),close(v$C_ri,b$C_ri,"central cri"))
      close(sum(v$C_ri),TOTAL_EXPECTED[[lab]],"printed central total",rel=ZERO,abs=ROUND_TOL)
    }
    ref_total <- NULL
    for (case in names(CASE_OMEGA)) {
      af <- read_pq(paste(case,"fi",lab,sep="_"))
      ac <- read_pq(paste(case,"cri",lab,sep="_"))
      need(af,c("grid_id","W_eff","f_i_full","SAR","SVR","p_d","k_used","p_l","fresh_fact","frac_sur_place"),case) # p_d replaces the unsaved p_d_safe.
      need(ac,c("grid_id","C_ri","C0i","di","f_i_full"),case)
      af <- align(af,a$grid_id,paste(case,lab,"fi")); ac <- align(ac,a$grid_id,paste(case,lab,"cri"))
      for (col in c("SAR","SVR","p_d","k_used","p_l","fresh_fact","frac_sur_place")) close(af[[col]],a[[col]],paste(case,col)) # close() compares NA masks too.
      close(ac$C0i,a$C0i,paste(case,"C0i")); close(ac$di,a$di,paste(case,"di"))
      omega <- CASE_OMEGA[[case]]
      for (arm in ARMS) {
        v <- calc(a,omega,PD_REF,arm)
        close(v$W_eff,af$W_eff,paste(case,"W"))
        close(v$fi_merge,af$f_i_full,paste(case,"fi"))
        close(v$fi,ac$f_i_full,paste(case,"step6 fi"))
        record(lab,paste(arm,case,"cellwise_cri"),close(v$C_ri,ac$C_ri,paste(case,"cri")))
        close(sum(v$C_ri),sum(ac$C_ri),paste(case,"total"),rel=TOTAL_REL,abs=ZERO)
        if (case == "ref") ref_total <- sum(ac$C_ri)
        if (is.null(ref_total) || ref_total <= ZERO) fail("Invalid Omega denominator")
        ref_calc <- sum(calc(a,OMEGA_REF,PD_REF,arm)$C_ri)
        ratio <- sum(v$C_ri)/ref_calc
        close(ratio,sum(ac$C_ri)/ref_total,paste(case,"unrounded ratio"),rel=TOTAL_REL,abs=ZERO)
        target <- EXPECTED$ratio[EXPECTED$label == lab & EXPECTED$omega == omega]
        record(lab,paste(arm,case,"published_ratio"),close(ratio,target,paste(case,"printed ratio"),rel=ZERO,abs=ROUND_TOL))
      }
    }
  }
  # GLOBAL GATE: no p_d other than 1 m calculated above, and both labels must pass.
  save_table(rbindlist(locks),"locks.tsv")
  write_new("PASS: all source/data/central/Omega locks, both labels and arms", "LOCK_PASSED")
  cat("LOCK PASSED: sensitivity calculations now permitted\n")

  rows <- list(); diags <- list()
  for (lab in LABELS) for (arm in ARMS) for (pd in PD) for (omega in OMEGA) {
    d <- data[[lab]]; v <- calc(d,omega,pd,arm)
    finite(v$C_ri,"scenario C_ri")
    tot <- sum(v$C_ri)
    # Count exact attained bounds, no arbitrary near-cap tolerance. Exclude zero stock from C cap.
    fcap <- !is.na(v$raw_f) & v$raw_f >= ONE
    cf <- v$C0i > ZERO & v$cpre >= v$C0i
    share <- function(mask) if(tot > ZERO) sum(v$C_ri[mask])/tot else NA_real_
    uncapped <- v$raw_f; uncapped[is.na(uncapped)] <- ZERO # Same NA handling as step6, for diagnostic only.
    rows[[length(rows)+1L]] <- data.table(label=lab,arm=arm,p_d_m=pd,omega_cm_yr=omega,
      C_ri_t=tot,C_ri_uncapped_t=sum(v$C0i*uncapped))
    # Saturation reports all cells and positive-SAR cells, plus stock-weighted share.
    diags[[length(diags)+1L]] <- data.table(label=lab,arm=arm,p_d_m=pd,omega_cm_yr=omega,
      n_cells=nrow(v),n_active=sum(d$SAR>ZERO),n_positive_stock=sum(d$C0i>ZERO),
      n_fi_upper=sum(fcap),n_fi_upper_active=sum(fcap & d$SAR>ZERO),share_Cri_fi_upper=share(fcap),
      n_fi_strict_clip=sum(v$raw_f>ONE,na.rm=TRUE),n_fi_lower=sum(v$raw_f<=ZERO,na.rm=TRUE),
      n_fi_NA_to_zero=sum(is.na(v$fi_merge)),n_Cri_at_stock=sum(cf),share_Cri_at_stock=share(cf),
      n_Cri_strict_clip=sum(v$cpre>v$C0i),
      n_frac_zero=sum(v$frac==ZERO),n_new_frac_zero=sum(v$frac==ZERO & d$frac_sur_place>ZERO),
      n_frac_missing_default=sum(is.na(d$V_retire_test)),
      C0i_frac_zero=sum(d$C0i[v$frac==ZERO]),
      max_fi_raw=if(any(!is.na(v$raw_f))) max(v$raw_f,na.rm=TRUE) else NA_real_)
  }
  res <- rbindlist(rows); diagnostics <- rbindlist(diags)
  den <- res[omega_cm_yr==OMEGA_REF,.(label,arm,p_d_m,denom_t=C_ri_t,denom_uncapped=C_ri_uncapped_t)]
  res <- merge(res,den,by=c("label","arm","p_d_m"),all.x=TRUE)
  res[, R:=fifelse(denom_t>ZERO,C_ri_t/denom_t,NA_real_)]
  res[, R_uncapped:=fifelse(denom_uncapped>ZERO,C_ri_uncapped_t/denom_uncapped,NA_real_)]
  baseline <- res[p_d_m==PD_REF,.(label,arm,omega_cm_yr,R_at_1m=R)]
  inter <- merge(res,baseline,by=c("label","arm","omega_cm_yr"))
  inter[, `:=`(I=R/R_at_1m-ONE, Delta_R=R-R_at_1m, E=R-ONE,
    Q_excess=fifelse(abs(R_at_1m-ONE)>ABS,(R-ONE)/(R_at_1m-ONE),NA_real_))]
  arms <- merge(res[arm=="frac_fixe",.(label,p_d_m,omega_cm_yr,R_fixed=R,C_fixed=C_ri_t)],
                res[arm=="frac_recalculee",.(label,p_d_m,omega_cm_yr,R_propagated=R,C_propagated=C_ri_t)],
                by=c("label","p_d_m","omega_cm_yr"))
  arms[, `:=`(B=R_propagated/R_fixed-ONE,Delta_R=R_propagated-R_fixed)]
  verdicts <- inter[, {
    trials <- .SD[omega_cm_yr!=OMEGA_REF & p_d_m!=PD_REF]
    primary <- .SD[omega_cm_yr==OMEGA_HIGH & p_d_m==min(PD)]
    valid <- all(is.finite(trials$I))
    max_i <- if(valid) max(abs(trials$I)) else NA_real_
    main_i <- primary$I
    list(transport=if(!valid) "NULL_DENOMINATOR" else if(max_i<MAT) "MAINTAINED_ON_GRID" else "CONDITIONAL_ON_DEPTH",
         max_abs_I=max_i, primary_I=main_i, primary_R=primary$R, primary_Q=primary$Q_excess,
         primary_transport=if(!is.finite(main_i)) "NULL_DENOMINATOR" else if(abs(main_i)<MAT) "MAINTAINED" else "FAILED_TRANSPORT",
         high_omega_at_30cm=if(!is.finite(primary$E)) "NULL_DENOMINATOR" else
           if(abs(primary$E)<MAT) "SECONDARY_LT_5PCT" else if(primary$E<=-MAT) "REVERSED" else "POSITIVE_GE_5PCT")
  },by=.(label,arm)]
  save_table(res,"results.tsv"); save_table(inter,"interaction.tsv")
  save_table(arms,"arms.tsv"); save_table(diagnostics,"diagnostics.tsv"); save_table(verdicts,"verdicts.tsv")
  save_table(inputs,"provenance.tsv")
  write_new(capture.output(sessionInfo()),"sessionInfo.txt")
  report <- c("# Omega x p_d: raw report (no manuscript edit)",
    paste("Slurm job:",Sys.getenv("SLURM_JOB_ID")),paste("Rule md5:",unname(tools::md5sum(file.path(T,"RULE.md")))),
    "All global locks PASS. Values below are conditional on fixed C0i over the upper metre and the steady profile.",
    "R denominator: same label, arm and excavation depth at reference omega. NA means a null comparison.",
    "Interaction I=R(p)/R(1m)-1; material at abs(I)>=0.05. Excess retention Q=(R(p)-1)/(R(1m)-1).",
    "Hopper contrast B=R_propagated/R_fixed-1; material at abs(B)>=0.05.",
    "Strict C cap counts are separate from exact equality at C0i; zero stock excluded from attained-cap counts.",
    capture.output(print(verdicts)),capture.output(print(res)),
    "Inspect diagnostics.tsv and arms.tsv before attribution. No stochastic inference or empirical parameter validation.")
  write_new(report,"RESULT.md")
  journal <- c("## Omega x p_d — execution record", "1. Question: Does the Omega effect depend on excavation depth?",
    "2. Preregistered rule: copy the archived run_claim/RULE.md; exact fingerprint and run below.",
    "3. Raw results:",report,"4. Verdict against original rule:",capture.output(print(verdicts)),
    "5. Provenance: provenance.tsv, bundle.md5 and results.md5; columns/rows identified by grid_id.",
    "6. Manuscript permissions: apply RULE.md outcome table. No wording is automatically authorized.")
  write_new(journal,"LOG_ENTRY.md")
  files <- list.files(OUT,pattern="[.](tsv|md|txt)$",full.names=TRUE)
  write_new(paste(unname(tools::md5sum(files)),basename(files),sep="  "),"results.md5")
  write_new("PASS: computation complete; copy the log record before interpreting", "COMPLETE")
  cat(paste(report,collapse="\n"),"\n")
}
tryCatch(main(),error=function(e) {
  message <- conditionMessage(e)
  if(grepl("token|secret|key|password|client_id|auth",message,ignore.case=TRUE)) message <- "Filtered diagnostic; inspect guards without exposing sensitive text"
  if(!file.exists(file.path(OUT,"locks.tsv")) && length(locks)) save_table(rbindlist(locks),"locks.tsv")
  text <- c("# Execution failed — stop interpretation",paste("Slurm job:",Sys.getenv("SLURM_JOB_ID")),
    message,"1. Question: Can the instrument reproduce the reference and complete the predefined test?",
    "2. Rule: RULE.md; any failed guard stops interpretation.","3. Raw result: diagnostic above and locks.tsv if present.",
    "4. Verdict: FAIL_INSTRUMENT_OR_INCOMPLETE; no sensitivity conclusion.",
    "5. Provenance: frozen inputs.tsv and bundle.md5; job identifier above.",
    "6. Manuscript: nothing. Copy this failure into the project log before any repair.")
  if(!file.exists(file.path(OUT,"FAILURE.md"))) write_new(text,"FAILURE.md")
  cat(paste(text,collapse="\n"),"\n")
  quit(status=1L)
})
