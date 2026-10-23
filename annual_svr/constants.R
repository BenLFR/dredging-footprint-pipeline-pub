# ────────────────────────────────────────────────────────────────────────────────
# SHARED CONSTANTS - Modular Step-5 pipeline
# 1 km Equal-Earth world grid, EPSG:6933
# ────────────────────────────────────────────────────────────────────────────────

# Global extent in metres (EPSG:6933) – consistent with Atwood and Step 3
WORLD_XMIN <- -17367530.45
WORLD_XMAX <-  17367530.45
WORLD_YMIN <- -7342699.72
WORLD_YMAX <-  7342699.72

# 1 km grid parameters
CELL_SIZE_M    <- 1000
CELL_AREA_M2   <- CELL_SIZE_M * CELL_SIZE_M
GRID_COLS      <- as.integer((WORLD_XMAX - WORLD_XMIN) / CELL_SIZE_M)   # 34735
GRID_ROWS      <- as.integer((WORLD_YMAX - WORLD_YMIN) / CELL_SIZE_M)   # 14685
NCOLS          <- GRID_COLS
NROWS          <- GRID_ROWS

# Tile buffers (SLURM multi-run test)
BUFFER_TEST_VALUES <- c(2000, 10000, 20000, 50000, 100000)
BUFFER_IDX <- as.integer(Sys.getenv("BUFFER_IDX", "4"))
if(is.na(BUFFER_IDX) || BUFFER_IDX < 1 || BUFFER_IDX > length(BUFFER_TEST_VALUES)) BUFFER_IDX <- 1
TILE_BUFFER_M <- BUFFER_TEST_VALUES[BUFFER_IDX]
cat(sprintf("🌍 [INFO] TILE_BUFFER_M set to %d m (BUFFER_IDX = %d)\n", TILE_BUFFER_M, BUFFER_IDX))

# Minimum arrow version (adjust to your software stack)
PARQUET_VERSION_MIN <- "14.0.0"

# Tile generation parameters
CELL_KM <- 1000  # Tile size in km
CRS_EQUIVALENT <- "EPSG:6933"  # Equal-Earth CRS

# Parameters for computing the effective lithology (thesis §2.7-2.8)
SURF_HORIZON <- 0.05          # 5 cm (surface horizon)
DEEP_HORIZON <- 0.05          # 5 cm max (deep horizon - CAPPED)
FACTOR_DEEP  <- 0.3           # Deep layer weighting factor

# --- Effective width of the scraping tool -----------------------------------
# `dredge_width_m` in ship_specs.yaml is the SUCTION PIPE DIAMETER
# (1.2-1.3 m; identical to suction_pipe_diameter_mm / 1000), not the width
# of the draghead that actually scrapes the seabed. Vlasblom, Designing Dredging
# Equipment, gives a draghead width / pipe diameter ratio
# of 3.0 (fig. 2-36; 2.8-3.06 over 450-1200 mm).
#
# The width enters the chain at ONE SINGLE point (W_v in step5), in a
# strictly linear way: SAR -> SVR -> f_i -> C_ri.
#
# Linearity domain, measured on fi_grid_20260328:
# the most dredged cell of the grid has f_i = 3.18 %, so the cap
# f_i <= 1 is only reached at a width multiplier of about 31.
# At 6.0 (option C): 0 saturated cells, deviation from linearity 0.00 %.
#
# The second cap, C_ri <= C0i, is STRUCTURALLY REDUNDANT:
# C_ri = pmin(C0i * f_i * di, C0i) with f_i <= 1 and di in {1.0; 0.272},
# so f_i * di <= 1 always. This pmin can never bind. Harmless,
# but it suggests a protection that is not one.
#
# WARNING: W_v is not
#     dredge_width_m * DRAGHEAD_WIDTH_FACTOR
# but
#     diam_m * n_pipes_eff * DRAGHEAD_WIDTH_FACTOR
# where n_pipes_eff = n_suction_pipes - (n_suction_pipes - 1) * PIPE_MONO_FRACTION.
# The factor below remains the draghead / pipe ratio of ONE head; the
# number of heads lives in ship_specs.yaml, NOT here.
#
# SET TO 1.0 TO RECOVER EXACTLY THE PUBLISHED BEHAVIOUR.
# Can be overridden without editing the code:  export DRAGHEAD_WIDTH_FACTOR=1.0
DRAGHEAD_WIDTH_FACTOR <- as.numeric(Sys.getenv("DRAGHEAD_WIDTH_FACTOR", "3.0"))
if (!is.finite(DRAGHEAD_WIDTH_FACTOR) || DRAGHEAD_WIDTH_FACTOR <= 0)
  stop("DRAGHEAD_WIDTH_FACTOR invalid: ", Sys.getenv("DRAGHEAD_WIDTH_FACTOR"))


# --- Cutting depth per pass ------------------------------------------------
# Depth of sediment removed by the draghead at each pass.
#
# 🔵 THIS VALUE IS A DECLARED TRUNCATION, NOT A MEASUREMENT. 1.0 m is
# the horizon to which Atwood et al. provide carbon data: beyond it,
# nothing can be accounted for. It is NOT the physical cutting
# depth, which may exceed it.
# The only published value close to it is 0.45 m
# (Miedema, Cristobal Colon case) -- but taken with its width of 4.26 it
# makes the mass balance impossible.
#
# The mass balance does not constrain width and depth separately,
# only their PRODUCT, and only FROM BELOW:
#        lambda x CUT_DEPTH_M >= 3.37 to 4.98
# A LOWER bound never refutes a deeper cut: the balance
# can impose a truncation, never exclude it.
# With one head the product is 3.0: the balance DOES NOT CLOSE.
# Decisive threshold: single-head time fraction <= 34 %.
# Can be overridden:  export CUT_DEPTH_M=0.45
CUT_DEPTH_M <- as.numeric(Sys.getenv("CUT_DEPTH_M", "1.0"))
if (!is.finite(CUT_DEPTH_M) || CUT_DEPTH_M <= 0)
  stop("CUT_DEPTH_M invalid: ", Sys.getenv("CUT_DEPTH_M"))

