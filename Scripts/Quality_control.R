source("~/Documents/Scripts/Ocelloscope/Growth_rate_calculations.R")

#Figure of distribution of R2 coefficients
r2hist <- ggplot(gr.res, aes(x = R2)) +
  geom_histogram(fill = "grey", colour = "black") +
  geom_vline(xintercept = 0.9, lty = "dashed", col = "red") +
  scale_y_continuous(expand = c(0,0))

#Distribution of number of tips
tipshist <- ggplot(gr.res, aes(x = tips.at.limit1)) +
  geom_histogram(fill = "grey", colour = "black") +
  scale_y_continuous(expand = c(0,0)) +
  xlab("No. tips in window")    

#Tips vs R2
tips.r2 <- ggplot(gr.res, aes(x = tips.at.limit1, y = R2)) +
  geom_point() +
  xlab("No. tips in window")     

finalplot <- plot_grid(r2hist, tipshist, tips.r2)
cowplot::save_plot(filename = rstudioapi::selectFile(caption = "Save diagnostic plot", label = "Save", filter = "PDF Files (*.pdf)", existing = FALSE), plot = finalplot)

summarize(gr.res, gr.ok = sum(R2 > 0.9, na.rm = TRUE), missing = sum(is.na(gr)), gr.bad = sum(R2 < 0.9, na.rm = TRUE))

