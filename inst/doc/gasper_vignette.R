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
knitr::opts_knit$set(global.par = TRUE)

## ---- eval=F------------------------------------------------------------------
#  head(SuiteSparseData, 15)

## ---- echo=FALSE--------------------------------------------------------------
SuiteSparseData_subset <- head(SuiteSparseData, 15)
df1 <-kableExtra::kbl(SuiteSparseData_subset, 
                    valign = 't',
                    format = "latex", booktabs = T,
                    caption="Overview of the first 15 matrices from the SuiteSparse Matrix Collection.")
kableExtra::kable_styling(df1, latex_options = c("striped","HOLD_position"))

## ---- eval=FALSE--------------------------------------------------------------
#  filtered_mat <- SuiteSparseData[SuiteSparseData$Kind == "Undirected Graph" &
#                   SuiteSparseData$Rows >= 100 & SuiteSparseData$Rows <= 150 &
#                   SuiteSparseData$Cols >= 100 & SuiteSparseData$Cols <= 150, ]
#  filtered_mat

## ---- echo=FALSE--------------------------------------------------------------
filtered_mat <- SuiteSparseData[SuiteSparseData$Kind == "Undirected Graph" & 
                 SuiteSparseData$Rows >= 100 & SuiteSparseData$Rows <= 150 &
                 SuiteSparseData$Cols >= 100 & SuiteSparseData$Cols <= 150, ]
df2 <-kableExtra::kbl(filtered_mat, 
                    valign = 't',
                    format = "latex", booktabs = T,
                    caption = "Subset of undirected matrices with 100 to 150 rows and columns.")
kableExtra::kable_styling(df2, latex_options = c("striped","HOLD_position"))

## ---- include=FALSE-----------------------------------------------------------
data(grid1)
grid1_1 <- grid1
grid1$info <- NULL

## ---- eval = F----------------------------------------------------------------
#  matrixname <- "grid1"
#  groupname <- "AG-Monien"
#  download_graph(matrixname, groupname)

## -----------------------------------------------------------------------------
attributes(grid1)

## -----------------------------------------------------------------------------
str(grid1$sA)

## -----------------------------------------------------------------------------
head(grid1$xy, 3)

## -----------------------------------------------------------------------------
grid1$dim

## ---- eval=FALSE--------------------------------------------------------------
#  list.files(grid1$temp)

## ---- eval = F----------------------------------------------------------------
#  file.show(paste(grid1$temp,"grid1",sep=""))

## ---- eval = F----------------------------------------------------------------
#  cat(readLines(paste(grid1$temp,"grid1",sep=""), n=14), sep = "\n")

## ---- eval=FALSE--------------------------------------------------------------
#  matrix_info <- get_graph_info(matrixname, groupname)
#  matrix_info

## ---- eval=FALSE--------------------------------------------------------------
#  downloaded_graph <- download_graph(matrixname, groupname, add_info = TRUE)
#  downloaded_graph$info

## ---- echo=FALSE, eval=T,fig.pos='H'------------------------------------------
data(grid1)
df1 <-kableExtra::kbl(grid1_1$info[[1]], 
                    valign = 't',
                    format = "latex", booktabs = T)
#kableExtra::kable_styling(df1, latex_options = c("striped","HOLD_position"))
df2<-kableExtra::kable(grid1_1$info[[2]], 
                 valign = 't',
                 format = "latex", booktabs = T)
#kableExtra::kable_styling(df2, latex_options = c("striped","HOLD_position"))
x_table <- knitr::kables(
  list(df1,df2), 
  format = "latex", 
  caption = "Matrix Information (left) and Matrix Properties (right).")
kableExtra::kable_styling(x_table,
                          latex_options = c("HOLD_position"))

## ---- fig.show='hold', fig.cap="Graph (left) and graph signal (right)."-------
f <- rnorm(nrow(grid1$sA))
plot_graph(grid1)
plot_signal(grid1, f, size = 2)

## ---- eval=FALSE, include=FALSE-----------------------------------------------
#  grid1$xy <- NULL
#  attributes(grid1)
#  plot_graph(grid1)

## -----------------------------------------------------------------------------
grid1$info$MatrixProp["Strongly Connect Components",]

## -----------------------------------------------------------------------------
W <- matrix(c(0, 1, 0,
              1, 0, 1,
              0, 1, 0), ncol=3)
laplacian_mat(W, "unnormalized")
laplacian_mat(W, "normalized")
laplacian_mat(W, "randomwalk")

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

## ---- echo=FALSE--------------------------------------------------------------
par(mgp = c(1.5, 0.5, 0), 
    tcl = 0.2,
    mar = .1 + c(2.5,2.5,0,0), 
    oma = c(0,0,0,0),
    cex.axis = 0.8,
    las = 1)

## ----fig.width=4,fig.height=2, fig.cap="Plot of the spectral graph filters on the spectrum of grid1 graph."----
plot_filter(lmax,b)

## -----------------------------------------------------------------------------
n <- nrow(L)
f <- randsignal(0.01, 3, A)
sigma <- 0.01
noise <- rnorm(n, sd = sigma)
tilde_f <- f + noise

## ---- fig.show='hold', fig.cap="Original graph signal (left) and noisy observations (right)."----
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

## ---- fig.width=4, fig.height=3, fig.cap="MSE risk and its SURE estimates as a function of the threshold parameter."----
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

## ----risk, echo=F-------------------------------------------------------------
df <-kableExtra::kbl(res, 
                    valign = 't',
                    format = "latex", booktabs = T, caption="Comparison of SNR performance between  uniform and dependent policies." )
kableExtra::kable_styling(df, latex_options = c("striped","HOLD_position"))
#knitr::kable(res, caption="Uniform vs Dependent" )

## ---- eval=FALSE--------------------------------------------------------------
#  hatf_oracle_u <- inverse_sgwt(wc_oracle_u,
#                                evalues, evectors, b)

## ---- eval=FALSE--------------------------------------------------------------
#  wc_oracle_u <- betathresh(wcn,
#                            thresh[opt_thresh_u$min[[1]]], 2)

## -----------------------------------------------------------------------------
J <- floor(log(lmax)/log(b)) + 2
LD_opt_thresh_u <- LD_SUREthresh(J=J, 
                             wcn=wcn, 
                             diagWWt=diagWWt, 
                             beta=2, 
                             sigma=sigma,
                             hatsigma=NA,
                             policy = "uniform",
                             keepSURE = FALSE)
hatf_LD_SURE_u <- synthesis(LD_opt_thresh_u$wcLDSURE, tf)
print(paste0("LD_SURE_u = ",round(SNR(f,hatf_LD_SURE_u),2),"dB"))

