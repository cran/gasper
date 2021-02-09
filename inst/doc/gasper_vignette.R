## ----setup, include = FALSE---------------------------------------------------
library(gasper)
if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}
set.seed(434343)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=3,
  fig.height=3, 
  fig.align="center"
)

## -----------------------------------------------------------------------------
graphname <- "grid1"
groupname <- "AG-Monien"
download_graph(graphname, groupname)
attributes(grid1)

## -----------------------------------------------------------------------------
str(grid1$sA)

## -----------------------------------------------------------------------------
head(grid1$xy, 3)

## -----------------------------------------------------------------------------
grid1$dim

## -----------------------------------------------------------------------------
cat(readLines(grid1$info, n=14), sep = "\n")

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

## -----------------------------------------------------------------------------
# wcf <- forward_sgwt(f, evalues, evectors, b=b)

## -----------------------------------------------------------------------------
diagWWt <- colSums(t(tf)^2)
thresh <- sort(abs(wcn))
opt_thresh <- SURE_MSEthresh(wcn, 
                           wcf, 
                           thresh, 
                           diagWWt, 
                           b=2, 
                           sigma, 
                           NA,
                           policy = "dependent")

## -----------------------------------------------------------------------------
plot(thresh, opt_thresh$res$MSE,
     type="l", xlab = "t", ylab = "risk")
lines(thresh, opt_thresh$res$SURE-n*sigma^2, col="red")
legend("topleft", legend=c("MSE", "SURE"),
       col=c("black", "red"),lty=1)

## -----------------------------------------------------------------------------
hatf_oracle <- synthesis(opt_thresh$wc[,opt_thresh$min[1]], tf)
hatf_SURE  <- synthesis(opt_thresh$wc[,opt_thresh$min[2]], tf)

SNR(f,tilde_f)
SNR(f,hatf_oracle)
SNR(f,hatf_SURE)

## -----------------------------------------------------------------------------
#hatf_oracle <- inverse_sgwt(opt_thresh$wc[,opt_thresh$min[1]], evalues, evectors, b)

