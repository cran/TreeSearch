## ----plot-concavities, fig.align="center", fig.width = 7, fig.height = 5, echo = FALSE----
e <- 0:15

IWTotal <- function (k) e / (e + k)
IW <- function (k) {
  iwTotal <- IWTotal(k)
  iwStep <- iwTotal[-1] - iwTotal[-length(e)]
  iwStep / max(iwStep)
}
IWPoints <- function (k, n) {
  points(IW(k) ~ e[-1], type = 'b', col = n, pch = 2, lty = 2)
}

plot(e[-1], e[-1], type = 'n', frame.plot = FALSE,
     xlab = 'Extra step',
     ylab = "Relative cost of this extra step",
     ylim = c(0, 1))
k <- c(3, 10, 30)
for (i in seq_along(k)) IWPoints(k[i], i)

PP <- function (a, col, n = 40) {
  char <- rep(0:1, c(a, n - a))
  si <- TreeSearch::StepInformation(char)
  points(I(si / max(si)) ~ I(seq_along(si)), type = 'b', col = col, pch = 1)
}
p <- c(3, 8, 20)
for (i in seq_along(p)) PP(p[i], i + 3)
legend('topright', col = 1:6, lty = c(rep(2, 3), rep(1, 3)),
       pch = c(rep(2, 3), rep(1, 3)),
       bty = 'n', cex = 0.9,
       c(paste0('IW, k = ', k),
         paste0('PP, ', p, ' | ', 40 - p)))



