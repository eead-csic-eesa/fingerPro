# Created by Ivan Lizaga and Borja Latorre 09/10/2020
#' Input sediment sources
#'
#' The function select and extract the source samples and run MixSIAR model.
#'
#' @param data Data frame containing source and mixtures data
#'
#' @export
#' 
MixSIAR_F_SD <- function(data, resid.err = F, process.err = T, Means = F, prior=1 , chain = "test" ) {
  library(tidyr)
  library(fingerPro)
run <- chain

if (Means == T) {
  
sources <-data[c(1:(nrow(data)-1)),c(2:ncol(data))]
x <- round((ncol(data)-3)/2 + 2)
mixtures <- inputSample(data[nrow(data),1:x])
colnames(sources)[1] <- "id"
# If it's here only apply the SD reduction when using means and SD 
# It should not be behind if you want to replicate the Consensus results 
# because there weren't done it with the reduced SD
n_SD_reduce <- (ncol(sources)-2)/2 + 2
sources[,c(n_SD_reduce:(ncol(sources)-1))] <- sources[,c(n_SD_reduce:(ncol(sources)-1))]/sqrt(sources[,ncol(sources)])

} else{
  ###########################################
  sources <- inputSource(data)
  mixtures <- inputSample(data)
  ###########################################
}

# n_SD_reduce <- (ncol(sources)-2)/2 + 2
# sources[,c(n_SD_reduce:(ncol(sources)-1))] <- sources[,c(n_SD_reduce:(ncol(sources)-1))]/sqrt(sources[,ncol(sources)])

source_names <- unique(data[2])
source_names<- source_names[-c(nrow(source_names)),]

sources[1] <-source_names

sources[sources ==0] <- 0.001


mixtures[1] <- NULL
write.csv(mixtures, "tmp-mix.csv", row.names = FALSE, quote = FALSE)

colnames(sources)[1] <- "sources"
n = (length(colnames(sources)) - 2 ) / 2
ids <- c(1)
for (i in seq(1, n))
{
  colnames(sources)[1+i] <- paste0("Mean", colnames(sources)[1+i])
  colnames(sources)[1+i+n] <- substring(colnames(sources)[1+i+n], 2)
  colnames(sources)[1+i+n] <- paste0("SD", colnames(sources)[1+i+n])
  ids <- c(ids, 1+i, 1+i+n)
}
ids <- c(ids, 1+2*n+1)
sources <- sources[,ids]

write.csv(sources, "tmp-sources.csv", row.names = FALSE, quote = FALSE)

##################################
# Adding discrimination data == 0
##################################
for (i in seq(1, n))
{
  sources[,1+i] <- 0
  sources[,1+i+n] <- 0
}
write.csv(sources, "tmp-disc.csv", row.names = FALSE, quote = FALSE)

file <- "tmp-mix.csv"
cols <- colnames(read.csv(file))
mix <- load_mix_data(
  filename=file, iso_names=cols,
  factors=NULL, fac_random=NULL,
  fac_nested=NULL, cont_effects=NULL
)

sources <- load_source_data(
  filename="tmp-sources.csv",
  ########################################################################
  # T or F when using or not concentration dependence data
  ########################################################################
  source_factors=NULL, conc_dep= F,
  data_type="means",	mix
)

discr.all <- load_discr_data(
  filename="tmp-disc.csv", mix
)

all.filename <- "tmp.jags"
all.resid.err <- resid.err
all.process.err <- process.err
write_JAGS_model(
  all.filename, all.resid.err,
  all.process.err, mix, sources
)

jags.m3all <- run_model(
  run=run, mix, sources, discr.all,
  all.filename, alpha.prior= prior,
  all.resid.err, all.process.err)

# test, very short, short, normal, long, very long, extreme

file.remove("tmp.jags")
file.remove("tmp-mix.csv")
file.remove("tmp-sources.csv")
file.remove("tmp-disc.csv")

getQuant <- function(x) quantile(x,probs=c(.025,.05,.25,.5,.75,.95,.975))

global_quants <- t(round(apply(jags.m3all$BUGSoutput$sims.list$p.global,2,getQuant),3))

global_quants[,4] # 50% median
# dir.create("1202")

mixsiar_1 <- jags.m3all$BUGSoutput$sims.list$p.global
# 
output_JAGS(jags.m3all, mix, sources)

result_MixSIAR<<-mixsiar_1
}
# output_JAGS(jags.m3all, mix, sources, output_options=list(summary_save = TRUE, summary_name = "summary_statistics", sup_post = FALSE, plot_post_save_pdf = FALSE, plot_post_name = "posterior_density", sup_pairs = FALSE, plot_pairs_save_pdf = FALSE, plot_pairs_name = "pairs_plot", sup_xy = TRUE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot", gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = FALSE, diag_name = "diagnostics", indiv_effect = FALSE, plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png = FALSE))
