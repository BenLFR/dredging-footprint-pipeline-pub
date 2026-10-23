library(sensitivity)
cat(as.character(packageVersion("sensitivity")), "\n")
print(ls(asNamespace("sensitivity"), pattern = "Salt|Eff"))
print(sensitivity:::estim.sobolSalt)
