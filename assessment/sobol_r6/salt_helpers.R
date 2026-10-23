lib <- dirname(find.package("sensitivity"))
db <- tools::Rd_db("sensitivity", lib.loc = lib)
tools::Rd2txt(db[["sobolSalt.Rd"]], out = "salt_help.txt", options = list(underline_titles = FALSE))
