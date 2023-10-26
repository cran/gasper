## ----setup, include = FALSE---------------------------------------------------
library(gasper)
if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}
set.seed(434343)
data(grid1)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=3,
  fig.height=3, 
  fig.align="center"
)

## ---- eval = T----------------------------------------------------------------
matrixname <- "grid1"
groupname <- "AG-Monien"
download_graph(matrixname, groupname)
attributes(grid1)

## -----------------------------------------------------------------------------
str(grid1$sA)

## -----------------------------------------------------------------------------
head(grid1$xy, 3)

## -----------------------------------------------------------------------------
grid1$dim

## -----------------------------------------------------------------------------
list.files(grid1$temp)

## ---- eval = FALSE------------------------------------------------------------
#  cat(readLines(paste(grid1$temp,"grid1",sep=""), n=14), sep = "\n")

## ---- eval=FALSE--------------------------------------------------------------
#  matrix_info <- get_graph_info(matrixname, groupname)
#  matrix_info

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  graph_info <- get_graph_info(matrixname, groupname)
#  knitr::kables(list(knitr::kable(graph_info[[2]], valign = 't'),
#                     knitr::kable(graph_info[[3]], valign = 't')))

## ---- fig.show='hold'---------------------------------------------------------
f <- rnorm(nrow(grid1$sA))
plot_graph(grid1)
plot_signal(grid1, f, size = 2)

## -----------------------------------------------------------------------------
A <- grid1$sA
L <- laplacian_mat(A)
val1 <- eigensort(L)
evalues <- val1$evalues
evectors <- val1$evectors
#- largest eigenvalue
lmax <- max(evalues)
#- parameter that controls the scale number
b <- 2
tf <- tight_frame(evalues, evectors, b=b)

## ----fig.width=5--------------------------------------------------------------
plot_filter(lmax,b)

## -----------------------------------------------------------------------------
n <- nrow(L)
f <- randsignal(0.01, 3, A)
sigma <- 0.01
noise <- rnorm(n, sd = sigma)
tilde_f <- f + noise

## ----fig.show='hold'----------------------------------------------------------
plot_signal(grid1, f, size = 2)
plot_signal(grid1, tilde_f, size = 2)

## -----------------------------------------------------------------------------
wcn <- analysis(tilde_f,tf)
wcf <- analysis(f,tf)

## ---- eval=FALSE--------------------------------------------------------------
#  wcf <- forward_sgwt(f, evalues, evectors, b=b)

## -----------------------------------------------------------------------------
diagWWt <- colSums(t(tf)^2)
thresh <- sort(abs(wcn))
opt_thresh_d <- SURE_MSEthresh(wcn, 
                           wcf, 
                           thresh, 
                           diagWWt, 
                           beta=2, 
                           sigma, 
                           NA,
                           policy = "dependent",
                           keepwc = TRUE)

opt_thresh_u <- SURE_MSEthresh(wcn, 
                           wcf, 
                           thresh, 
                           diagWWt, 
                           beta=2, 
                           sigma, 
                           NA,
                           policy = "uniform",
                           keepwc = TRUE)

## ---- fig.width=5, fig.height=4-----------------------------------------------
plot(thresh, opt_thresh_u$res$MSE,
     type="l", xlab = "t", ylab = "risk", log="x")
lines(thresh, opt_thresh_u$res$SURE-n*sigma^2, col="red")
lines(thresh, opt_thresh_d$res$MSE, lty=2)
lines(thresh, opt_thresh_d$res$SURE-n*sigma^2, col="red", lty=2)
legend("topleft", legend=c("MSE_u", "SURE_u",
                           "MSE_d", "SURE_d"),
       col=rep(c("black", "red"), 2), 
       lty=c(1,1,2,2), cex = 1)

## -----------------------------------------------------------------------------
wc_oracle_u <- opt_thresh_u$wc[, opt_thresh_u$min["xminMSE"]]
wc_oracle_d <- opt_thresh_d$wc[, opt_thresh_d$min["xminMSE"]]
wc_SURE_u <- opt_thresh_u$wc[, opt_thresh_u$min["xminSURE"]]
wc_SURE_d <- opt_thresh_d$wc[, opt_thresh_d$min["xminSURE"]]

hatf_oracle_u <- synthesis(wc_oracle_u, tf)
hatf_oracle_d <- synthesis(wc_oracle_d, tf)
hatf_SURE_u  <- synthesis(wc_SURE_u, tf)
hatf_SURE_d  <- synthesis(wc_SURE_d, tf)

res <- data.frame("Input_SNR"=round(SNR(f,tilde_f),2),
                  "MSE_u"=round(SNR(f,hatf_oracle_u),2),
                  "SURE_u"=round(SNR(f,hatf_SURE_u),2),
                  "MSE_d"=round(SNR(f,hatf_oracle_d),2),
                  "SURE_d"=round(SNR(f,hatf_SURE_d),2))

## ---- echo=F------------------------------------------------------------------
knitr::kable(res, caption="Uniform vs Dependent" )

## ---- eval=FALSE--------------------------------------------------------------
#  hatf_oracle_u <- inverse_sgwt(wc_oracle_u,
#                                evalues, evectors, b)

## ---- eval=FALSE--------------------------------------------------------------
#  wc_oracle_u <- betathresh(wcn,
#                            thresh[opt_thresh_u$min[[1]]], 2)

## -----------------------------------------------------------------------------
J <- floor(log(lmax)/log(b)) + 2
LD_opt_thresh_d <- LD_SUREthresh(J=J, 
                             wcn=wcn, 
                             diagWWt=diagWWt, 
                             beta=2, 
                             sigma=sigma,
                             hatsigma=NA,
                             policy = "uniform",
                             keepSURE = FALSE)
hatf_LD_SURE_d <- synthesis(LD_opt_thresh_d$wcLDSURE, tf)
print(paste0("LD_SURE_u = ",round(SNR(f,hatf_LD_SURE_d),2),"dB"))

