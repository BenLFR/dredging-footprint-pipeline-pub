suppressMessages({library(data.table); library(arrow); library(terra)})
R <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_sala2021_r6_flux10")
carbon_dir <- path.expand(file.path(Sys.getenv("DFP_DATA"), "configuration/atwood_carbon_full"))
r <- rast(file.path(carbon_dir, "Mean carbon_stock.tif")); print(r)
for (lab in c("smooth","flag")) {
  f <- list.files(file.path(R,"cri_pr"), paste0("^cri_final_",lab,"_.*[.]parquet$"), full.names=TRUE)
  d <- as.data.table(read_parquet(f))[C0i == 0 & f_i_full > 0]
  pts <- project(vect(as.data.frame(d[, .(x, y)]), geom=c("x","y"), crs="EPSG:6933"), crs(r))
  v <- terra::extract(r, pts, ID=FALSE)[[1]]
  cat(lab, ": cells C0i==0 & f>0 =", nrow(d), " ; raster NA =", sum(is.na(v)), " ; raster 0 =", sum(v == 0, na.rm=TRUE),
      " ; other =", sum(!is.na(v) & v != 0), "\n")
}
