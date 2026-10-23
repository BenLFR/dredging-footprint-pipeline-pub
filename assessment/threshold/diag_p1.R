# diag_p1.R -- Diagnostic of the failed lock P1. OUTSIDE VERDICT: reads no C_ri,
# touches no arm. Question: where does the list of 71 tiles come from, and what does it lose?
#   (1) rows of each sar of the 71 tiles, flag and smooth, in _r2 (which tiles are empty?)
#   (2) tiles carrying flag / smooth / flag-or-smooth positions, with and without the 50 km buffer
#   (3) flag positions that none of the 71 enlarged tiles covers (positions never processed by step5)
# Usage: Rscript diag_p1.R
suppressMessages({ library(data.table); library(sf); library(arrow) })
R2    <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_sala2021_r2")
IN    <- file.path(R2, "AIS_with_lithology_20260922_221357.rds")
TILES <- path.expand(file.path(Sys.getenv("DFP_DATA"), "output_V6/tiles_1000km.gpkg"))
T71   <- c(104,105,118,125,153,154,169,185,186,189,197,205,206,220,232,238,239,241,255,262,263,270,272,273,276,
           289,290,291,296,297,298,305,306,307,308,323,327,331,337,338,339,340,341,344,345,357,366,371,372,373,374,
           379,380,401,402,403,404,405,406,415,437,438,439,440,441,451,452,473,475,511,515)  # chain_step4_ocim.sh
ec <- new.env(); invisible(capture.output(sys.source(path.expand(file.path(Sys.getenv("DFP_REPO"), "pipeline/constants.R")), envir = ec)))
cat("==== diag_p1.R ", format(Sys.time()), " ; s2 active:", sf_use_s2(), "\n")

# (1) sar of _r2
for (lab in c("flag", "smooth")) {
  f <- list.files(file.path(R2, paste0("step5_", lab)), "^sar_\\d+\\.parquet$", full.names = TRUE)
  id <- as.integer(sub("^sar_(\\d+)\\.parquet$", "\\1", basename(f)))
  n <- vapply(f, function(x) nrow(read_parquet(x)), 0L)
  cat(sprintf("\n(1) %s: %d sar ; = the 71: %s ; empty (0 rows): %s\n", lab, length(f), setequal(id, T71),
              paste(sort(id[n == 0]), collapse = ",")))
}

# (2) tiles carrying positions
a  <- readRDS(IN)
tl <- st_read(TILES, quiet = TRUE); ids <- tl$tile_id
cat(sprintf("\n    CRS of the tiles: %s\n", st_crs(tl)$input))
tlb <- st_buffer(st_geometry(tl), ec$TILE_BUFFER_M)
ok  <- !is.na(a[["Lon"]]) & !is.na(a[["Lat"]])
fl  <- a[["Dragage_flag"]] == 1 & ok
sm  <- a[["behavior_smooth"]] %in% "dredging" & ok
pts_of <- function(sel) { p <- unique(data.table(x = a[["Lon"]][sel], y = a[["Lat"]][sel]))
  list(p = p, sf = st_transform(st_as_sf(p, coords = c("x", "y"), crs = 4326), st_crs(tl))) }
tiles <- function(P, geom) sort(unique(ids[unlist(st_intersects(P$sf, geom))]))
P <- list(flag = pts_of(fl), smooth = pts_of(sm), union = pts_of(fl | sm), all_pts = pts_of(ok))
cat("\n(2) tiles carrying at least one position     without buffer | with 50 km buffer\n")
for (nm in names(P)) {
  s0 <- tiles(P[[nm]], st_geometry(tl)); s1 <- tiles(P[[nm]], tlb)
  cat(sprintf("    %-7s (%8d unique points): %3d (=71 %s) | %3d (=71 %s)\n", nm, nrow(P[[nm]]$p),
              length(s0), setequal(s0, T71), length(s1), setequal(s1, T71)))
  cat(sprintf("        without buffer: outside list %s ; missing %s\n", paste(setdiff(s0, T71), collapse = ","),
              paste(setdiff(T71, s0), collapse = ",")))
  cat(sprintf("        with buffer: outside list %s ; missing %s\n", paste(setdiff(s1, T71), collapse = ","),
              paste(setdiff(T71, s1), collapse = ",")))
}

# (3) flag positions that none of the 71 enlarged tiles covers
cov <- lengths(st_intersects(P$flag$sf, tlb[ids %in% T71])) > 0
key <- paste(a[["Lon"]], a[["Lat"]])
outside <- key %in% paste(P$flag$p$x[!cov], P$flag$p$y[!cov]) & fl
cat(sprintf("\n(3) flag positions outside the 71 enlarged tiles: %d of %d (%.4f %%) ; unique points %d\n",
            sum(outside), sum(fl), 100 * sum(outside) / sum(fl), sum(!cov)))
if (sum(!cov) > 0) {
  t_outside <- ids[unlist(st_intersects(P$flag$sf[!cov, ], st_geometry(tl)))]
  print(table(tile_without_buffer = t_outside, useNA = "ifany"))
  cat(sprintf("    extent: Lon %.3f to %.3f ; Lat %.3f to %.3f\n", min(P$flag$p$x[!cov]), max(P$flag$p$x[!cov]),
              min(P$flag$p$y[!cov]), max(P$flag$p$y[!cov])))
}
# same question for smooth
covs <- lengths(st_intersects(P$smooth$sf, tlb[ids %in% T71])) > 0
cat(sprintf("    smooth: unique points outside the 71 enlarged tiles: %d of %d\n", sum(!covs), length(covs)))
cat("END diag_p1\n")
