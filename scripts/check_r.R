# Parse source without loading packages or executing scientific calculations.
root <- commandArgs(TRUE)[1]
files <- unlist(lapply(c("pipeline", "config", "ocim", "annual_svr", "assessment", "scripts"), function(d)
  list.files(file.path(root, d), pattern="[.][Rr]$", recursive=TRUE, full.names=TRUE)))
for (f in files) parse(f, keep.source=FALSE, encoding="UTF-8")
cat("Syntax: R", length(files), "files parsed\n")
