
library(coloc)

args <- commandArgs(trailingOnly = TRUE)
in_csv <- args[1]
out_csv <- args[2]
type1 <- args[3]
type2 <- args[4]
N1 <- as.numeric(args[5])
N2 <- as.numeric(args[6])
s2 <- as.numeric(args[7])

df <- read.csv(in_csv)

d1 <- list(
  beta = df$beta_eqtl,
  varbeta = df$varbeta_eqtl,
  snp = df$snp,
  type = type1,
  N = N1
)

d2 <- list(
  beta = df$beta_gwas,
  varbeta = df$varbeta_gwas,
  snp = df$snp,
  type = type2,
  N = N2
)
if (!is.na(s2) && s2 > 0) {
  d2$s <- s2
}

res <- coloc.abf(dataset1 = d1, dataset2 = d2)
summary_res <- as.data.frame(t(res$summary))
write.csv(summary_res, out_csv, row.names = FALSE)
