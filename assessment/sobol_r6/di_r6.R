suppressMessages({library(data.table); library(arrow)})
R <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_sala2021_r6_flux10")
for (lab in c("smooth","flag")) {
  f <- list.files(file.path(R,"cri_pr"), paste0("^cri_final_",lab,"_.*[.]parquet$"), full.names=TRUE)
  d <- as.data.table(read_parquet(f))
  cat("==", lab, basename(f), nrow(d), "rows\n")
  print(d[, .N, by=di])
  x <- d[!is.na(C_ri) & !is.na(f_i_full)]
  r <- x[, C_ri - pmin(C0i*pmin(pmax(f_i_full,0),1)*di, C0i)]
  cat("max |C_ri - C0i f di| =", max(abs(r)), " ; sum C_ri =", sum(x$C_ri), "\n")
  cat("sum C0i f (without di) =", x[, sum(C0i*pmin(pmax(f_i_full,0),1))], "\n")
  cat("C0i: quantiles", quantile(d$C0i, c(0,.5,.99,1), na.rm=TRUE), "\n")
  cat("columns:", paste(names(d), collapse=" "), "\n")
}
