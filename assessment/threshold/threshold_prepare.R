#!/usr/bin/env Rscript
# threshold_prepare.R -- Classification threshold +/- 0.005 (Label B): preparation of the three arms.
# Rule written beforehand. Model: the earlier Label B preparation (p72_prepare.R, not published).
#
# Recomputes Dragage_flag as the installed step3 does:
#   flag = drag_score >= threshold (drag_score NA -> 0), then smooth_binary_runs(flag, Speed, thr = q75(n_min), 3.5)
#   per Seg_id. The smoothing function is PARSED from the installed step3, not copied (as in p72).
# Arms: t0 (stored threshold), tp (t0 + 0.005), tm (t0 - 0.005).
#
# Locks, in order; the first that fails stops the script (stop), nothing is written for step5:
#   counts of the input (8 235 295 / 1 837 000 / 1 948 150 / 1 464 644, p72 locks), run threshold = 6,
#   stored threshold = 0.0679029990731113 to 1e-12 (read in the original gridsearch, as in p72);
#   P0  at t0, the recomputed flag = the stored flag, position by position. File order first;
#       otherwise ONE alternative declared beforehand: sort by (Seg_id, Timestamp). The order that passes is used for the arms t+/-.
#   P1  REVISED before any arm was run (rule written beforehand). The 1st version ("the
#       tiles carrying >= 1 flag position = the 71") FAILED (66 tiles): the list of 71 comes from the
#       original step5 outputs, not from a rule (chain_step4_ocim.sh). Version in force:
#       - the three arms run on the FIXED list of 71, as in production;
#       - P1a (lock): at t0, no position Dragage_flag == 1 (Lon/Lat not NA) outside the union of the 71 tiles
#         enlarged by TILE_BUFFER_M (constants.R: 50 000 m, BUFFER_IDX not set -> 4; worker
#         step5_tile_worker_batch.R);
#       - P1b (reported, does not block): same count at t+ and t-; > 0 -> the Delta of this arm is declared incomplete;
#       - the guard "none of the 71 tiles emptied" is removed: production already runs 7 empty flag tiles.
#       Tiles carrying positions according to our geometry are printed, outside verdict.
#
# Outputs (BASE/<arm>/): AIS_with_lithology_seuil_<arm>.rds (t0: symbolic link to the _r2 input, not rewritten),
#   tuiles_<arm>.txt (list for --array); BASE/compteurs.csv, BASE/locks.csv, BASE/PREPARATION_PASS.
# Usage: Rscript threshold_prepare.R <BASE>
suppressMessages({ library(data.table); library(sf) })
BASE  <- commandArgs(TRUE)[1]
R2    <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_sala2021_r2")
IN    <- file.path(R2, "AIS_with_lithology_20260922_221357.rds")          # input of the step5 of _r2
GS    <- path.expand(file.path(Sys.getenv("DFP_DATA"), "output_V6_rejeu/dragage_gridsearch_results_V6_20260810_093400.rds"))  # as in p72_prepare.R
S3    <- path.expand(file.path(Sys.getenv("DFP_REPO"), "pipeline/step3_merge_final.R"))    # installed step3
TILES <- path.expand(file.path(Sys.getenv("DFP_DATA"), "output_V6/tiles_1000km.gpkg"))             # step5_batch.sbatch (PIPE_TILES)
T71   <- c(104,105,118,125,153,154,169,185,186,189,197,205,206,220,232,238,239,241,255,262,263,270,272,273,276,
           289,290,291,296,297,298,305,306,307,308,323,327,331,337,338,339,340,341,344,345,357,366,371,372,373,374,
           379,380,401,402,403,404,405,406,415,437,438,439,440,441,451,452,473,475,511,515)  # chain_step4_ocim.sh
T0_EXP <- 0.0679029990731113     # p72 lock, tolerance 1e-12 as in p72_prepare.R
DELTA  <- 0.005                  # design choice, no source
SPD_MAX <- 3.5                   # installed step3, line "spd_max = 3.5"
EXP <- c(n = 8235295, flag = 1837000, above = 1948150, smooth = 1464644)  # p72 locks
CONST <- path.expand(file.path(Sys.getenv("DFP_REPO"), "pipeline/constants.R"))           # source of the tile buffer
BUFFER_EXP <- 50000   # m; constants.R (BUFFER_IDX = 4 by default, not set by step5_batch)

dir.create(BASE, showWarnings = FALSE, recursive = TRUE)
locks <- data.table(check = character(), observed = character(), pass = logical())
lock <- function(name, value, ok) {
  locks <<- rbind(locks, data.table(check = name, observed = paste(value, collapse = ";"), pass = isTRUE(ok)))
  fwrite(locks, file.path(BASE, "locks.csv"))
  cat(sprintf("  %-5s %-34s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, paste(value, collapse = ";")))
  if (!isTRUE(ok)) stop("LOCK FAILED: ", name, " -- stop, nothing is written for step5")
}
cat("==== threshold_prepare.R ", format(Sys.time()), " ====\n")
lock("no earlier PREPARATION_PASS", BASE, !file.exists(file.path(BASE, "PREPARATION_PASS")))
lock("BUFFER_IDX not set in the environment", Sys.getenv("BUFFER_IDX"), !nzchar(Sys.getenv("BUFFER_IDX")))
ec <- new.env(); invisible(capture.output(sys.source(CONST, envir = ec)))
lock("tile buffer (constants.R)", ec$TILE_BUFFER_M, identical(as.numeric(ec$TILE_BUFFER_M), BUFFER_EXP))

# ---- 1. smoothing: parsed from the installed step3 (as in p72_prepare.R) --------------------------------
cat(sprintf("  installed step3: %s  md5 %s (outside lock ; expected md5 39617684 = sha256 857ff85a)\n",
            S3, unname(tools::md5sum(S3))))
src <- readLines(S3, warn = FALSE)
j <- grep("^        smooth_binary_runs <- function", src)
lock("smoother found in step3", length(j), length(j) == 1L)
e <- new.env(parent = baseenv()); done <- FALSE
for (k in j:length(src)) {
  ex <- tryCatch(parse(text = paste(src[j:k], collapse = "\n")), error = function(err) NULL)
  if (!is.null(ex)) { eval(ex, envir = e); done <- TRUE; break }
}
lock("smoother parsed", done, done)
lock("step3 : spd_max = 3.5", sum(grepl("spd_max = 3.5", src, fixed = TRUE)), any(grepl("spd_max = 3.5", src, fixed = TRUE)))
lock("step3: q75 of n_min", sum(grepl("quantile(ais$n_min, 0.75, na.rm = TRUE)", src, fixed = TRUE)),
     any(grepl("quantile(ais$n_min, 0.75, na.rm = TRUE)", src, fixed = TRUE)))

# ---- 2. input and counts ------------------------------------------------------------------------------
t0 <- readRDS(GS)$best_weights$optimal_threshold
lock("stored threshold", sprintf("%.16g", t0), abs(t0 - T0_EXP) < 1e-12)
cat(sprintf("  input md5: %s\n", unname(tools::md5sum(IN))))
a <- readRDS(IN)
gi <- attr(a, "gmm_info")
cat(sprintf("  class of the input: %s ; gmm_info present: %s\n", paste(class(a), collapse = ","), !is.null(gi)))
req <- c("Seg_id", "Timestamp", "Speed", "drag_score", "Dragage_flag", "n_min", "behavior_smooth", "Lon", "Lat")
lock("columns", setdiff(req, names(a)), all(req %in% names(a)))
orig <- a[["Dragage_flag"]]
sc <- a[["drag_score"]]
cat(sprintf("  drag_score NA: %d (set to 0, as in step3)\n", sum(is.na(sc))))
sc[is.na(sc)] <- 0
got <- c(n = length(orig), flag = sum(orig == 1, na.rm = TRUE), above = sum(sc >= t0),
         smooth = sum(a[["behavior_smooth"]] == "dredging", na.rm = TRUE))
for (nm in names(EXP)) lock(paste("count", nm), got[[nm]], got[[nm]] == EXP[[nm]])
lock("Dragage_flag without NA", sum(is.na(orig)), !anyNA(orig))
thr <- as.numeric(quantile(a[["n_min"]], 0.75, na.rm = TRUE))
lock("run threshold = 6", thr, thr == 6)

# ---- 3. P0: reproduction of the stored flag ---------------------------------------------------------------
smooth_flags <- function(t, ord) {
  d <- data.table(seg = a[["Seg_id"]][ord], spd = a[["Speed"]][ord], f = as.integer(sc[ord] >= t))
  d[, f := e$smooth_binary_runs(f, spd, thr = thr, spd_max = SPD_MAX), by = seg]
  out <- integer(length(ord)); out[ord] <- d$f; out
}
orders <- list(file_order = seq_along(orig),
               seg_time = order(a[["Seg_id"]], a[["Timestamp"]], seq_along(orig)))   # declared alternative
ord <- NULL
for (o in names(orders)) {
  f0 <- smooth_flags(t0, orders[[o]])
  nd <- sum(f0 != orig)
  cat(sprintf("  P0 order %-9s: %d differing positions (flag recomputed %d)\n", o, nd, sum(f0)))
  if (nd == 0L) { ord <- orders[[o]]; order_P0 <- o; break }
}
lock("P0 reproduction of the flag", if (is.null(ord)) "no order" else order_P0, !is.null(ord))

# ---- 4. flags of the three arms ---------------------------------------------------------------------------
thresholds <- c(t0 = t0, tp = t0 + DELTA, tm = t0 - DELTA)
flags <- list(t0 = f0)
for (b in c("tp", "tm")) flags[[b]] <- smooth_flags(thresholds[[b]], ord)
counts <- rbindlist(lapply(names(thresholds), function(b) data.table(
  bras = b, seuil = sprintf("%.16g", thresholds[[b]]), au_dessus = sum(sc >= thresholds[[b]]),
  flag = sum(flags[[b]]), flag_rel_t0_pct = 100 * (sum(flags[[b]]) / sum(f0) - 1),
  au_dessus_rel_t0_pct = 100 * (sum(sc >= thresholds[[b]]) / sum(sc >= t0) - 1),
  changees = sum(flags[[b]] != f0))))
print(counts)

# ---- 5. P1: tiles -----------------------------------------------------------------------------------------
tl <- st_read(TILES, quiet = TRUE)
cat(sprintf("  tiles: %d rows ; columns: %s\n", nrow(tl), paste(setdiff(names(tl), attr(tl, "sf_column")), collapse = ",")))
lock("tile_id column in the gpkg (worker)", "tile_id", "tile_id" %in% names(tl))
lock("tile_id unique", anyDuplicated(tl$tile_id), anyDuplicated(tl$tile_id) == 0L)
ids <- tl$tile_id
tlb <- st_buffer(st_geometry(tl), ec$TILE_BUFFER_M)             # worker (step5_tile_worker_batch.R)
ok_xy <- !is.na(a[["Lon"]]) & !is.na(a[["Lat"]])               # worker (step5_tile_worker_batch.R)
cat(sprintf("  positions without Lon/Lat: %d (discarded, as in the worker)\n", sum(!ok_xy)))
tiles_of <- function(fl, geom) {
  p <- unique(data.table(x = a[["Lon"]][fl == 1 & ok_xy], y = a[["Lat"]][fl == 1 & ok_xy]))
  pts <- st_transform(st_as_sf(p, coords = c("x", "y"), crs = 4326), st_crs(tl))
  sort(unique(ids[unlist(st_intersects(pts, geom))]))
}
# revised P1: coverage by the fixed list of 71
lock("the 71 tiles present in the gpkg", sum(ids %in% T71), sum(ids %in% T71) == 71L)
u71 <- st_union(tlb[ids %in% T71])                               # union of the 71 enlarged tiles
hors_71 <- function(fl) {                                        # positions classified outside the 71 enlarged tiles
  sel <- fl == 1 & ok_xy
  p <- unique(data.table(x = a[["Lon"]][sel], y = a[["Lat"]][sel]))
  pts <- st_transform(st_as_sf(p, coords = c("x", "y"), crs = 4326), st_crs(tl))
  out <- lengths(st_intersects(pts, u71)) == 0
  if (!any(out)) return(0L)
  sum(paste(a[["Lon"]][sel], a[["Lat"]][sel]) %in% paste(p$x[out], p$y[out]))
}
h <- vapply(flags, hors_71, 0L)
lock("P1a t0: no flag position outside the 71", h[["t0"]], h[["t0"]] == 0L)
for (b in c("tp", "tm")) {
  ch <- sum(flags[[b]] != f0)
  cat(sprintf("  P1b %s (reported): %d positions classified outside the 71 enlarged tiles (%.4f %% of the %d that change)%s\n",
              b, h[[b]], 100 * h[[b]] / max(ch, 1L), ch, if (h[[b]] > 0L) " -> Delta of this arm INCOMPLETE" else ""))
}
# outside verdict: tiles carrying positions according to our geometry (differs from the worker on tile 263, diag_p1)
uniq_tiles <- lapply(flags, tiles_of, geom = tlb)
cat(sprintf("  (outside verdict) t0 without buffer: %d tiles\n", length(tiles_of(flags$t0, st_geometry(tl)))))
for (b in names(flags))
  cat(sprintf("  (outside verdict) %s: %d tiles with buffer ; of the 71 without position: %s ; outside list with positions: %s\n",
              b, length(uniq_tiles[[b]]), paste(setdiff(T71, uniq_tiles[[b]]), collapse = ","), paste(setdiff(uniq_tiles[[b]], T71), collapse = ",")))
counts[, hors_71 := vapply(bras, function(b) h[[b]], 0L)]
counts[, tiles_geom := vapply(bras, function(b) length(uniq_tiles[[b]]), 1L)]
counts[, tiles_not_listed := vapply(bras, function(b) paste(setdiff(uniq_tiles[[b]], T71), collapse = " "), "")]
counts[, order_P0 := order_P0]
fwrite(counts, file.path(BASE, "compteurs.csv"))

# ---- 6. writing the inputs of step5 -----------------------------------------------------------------------
conv <- function(f) if (is.integer(orig)) as.integer(f) else if (is.double(orig)) as.double(f) else
  if (is.logical(orig)) as.logical(f) else stop("unexpected Dragage_flag type: ", class(orig))
for (b in names(flags)) {
  d <- file.path(BASE, b); dir.create(d, showWarnings = FALSE)
  writeLines(paste(T71, collapse = ","), file.path(d, sprintf("tuiles_%s.txt", b)))   # fixed list (revised P1)
  f <- file.path(d, sprintf("AIS_with_lithology_seuil_%s.rds", b))
  if (b == "t0") {
    lock("link t0 created", f, file.symlink(IN, f))
  } else {
    a[["Dragage_flag"]] <- conv(flags[[b]])
    lock(sprintf("%s: gmm_info kept", b), !is.null(attr(a, "gmm_info")), identical(attr(a, "gmm_info"), gi))
    lock(sprintf("%s: count written", b), sum(a[["Dragage_flag"]] == 1), sum(a[["Dragage_flag"]] == 1) == sum(flags[[b]]))
    saveRDS(a, f)
    cat(sprintf("  %s written: %s  md5 %s\n", b, f, unname(tools::md5sum(f))))
  }
}
writeLines(c("PREPARATION_PASS", format(Sys.time()), paste("P0 order:", order_P0), paste("buffer:", ec$TILE_BUFFER_M, "m")),
           file.path(BASE, "PREPARATION_PASS"))
cat("PREPARATION_PASS\n")
