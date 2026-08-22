log10p <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)

# LocusZoom style view: -log10(P) vs position, coloured by r^2 to the lead SNP.

locuszoom_plot <- function(ss, R, title = NULL) {
  lead <- which.max(abs(ss$Z))
  d <- data.table(POS = ss$POS, logp = log10p(ss$Z), r2 = R[lead, ]^2)
  ggplot(d, aes(POS / 1e6, logp, colour = r2)) +
    geom_point(size = 1.6) +
    geom_point(data = d[lead], colour = "purple", shape = 18, size = 4) +
    scale_colour_gradientn(colours = c("navy", "skyblue", "green", "orange", "red"),
                           limits = c(0, 1), name = expression(r^2)) +
    labs(x = "Position (Mb)", y = expression(-log[10](P)), title = title) +
    theme_bw()
}

# PIP track: posterior inclusion probability vs position, coloured by credible
# set.
pip_track <- function(fit, ss, title = NULL) {
  pip  <- susie_get_pip(fit)
  csid <- rep("none", length(pip))
  cs   <- fit$sets$cs
  if (length(cs)) for (i in seq_along(cs)) csid[cs[[i]]] <- paste0("CS", i)
  d <- data.table(POS = ss$POS, PIP = pip, CS = factor(csid))
  ggplot(d, aes(POS / 1e6, PIP, colour = CS)) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c(none = "grey75", CS1 = "firebrick",
                                   CS2 = "steelblue", CS3 = "darkgreen")) +
    ylim(0, 1) +
    labs(x = "Position (Mb)", y = "PIP", title = title) +
    theme_bw()
}