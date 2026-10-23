suppressPackageStartupMessages(library(arrow))
x <- read.delim(file.path(Sys.getenv("DFP_REPO"), "assessment/omega_pd/inputs.tsv"), colClasses="character")
fi <- c("grid_id","SAR","SVR","p_d","p_l","k_used","fresh_fact","W_eff","f_i_full","f_i_full_avant_pr","frac_sur_place")
cr <- c("grid_id","C0i","C0i_mean_carbon_stock","di","f_i_full","C_ri")
for (i in which(grepl("[.]parquet$", x$path))) {
  nm <- names(open_dataset(x$path[i])$schema)
  req <- if (grepl("_fi_", x$id[i])) fi else cr
  miss <- setdiff(req, nm)
  cat(x$id[i], "| missing:", if (length(miss)) paste(miss, collapse=",") else "none", "| p_d_safe present:", "p_d_safe" %in% nm, "\n")
}
