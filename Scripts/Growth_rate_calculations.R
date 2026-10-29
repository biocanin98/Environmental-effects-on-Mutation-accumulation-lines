source("oCelloscope_functions.R")

#Choose file automatically
datafile <- file.choose()

#Load data file
testdata <- load.ocelloscope.data2(datafile, sep = ",", dec = ".", repetitions = 60, features = c("EstimatedVolume", "HyphaeArea", "NumberOfTips", "TipsPerArea", "TimestampInSeconds"))

#Convert to longdata
longdata <- convert.ocelloscope.data.long(testdata)

#Calculate growth rate
gr.res <- ocelloscope.growth.rate(longdata)

# Plot the curves
ggplot(longdata, aes(x = Time, y = area, group = sample)) +
  geom_line() +
  geom_vline(xintercept = c(0, 15000), lty = "dashed", col = "red") +
  ylab("Area")

# Isolate windowed area
  ggplot(
    dplyr::filter(longdata, Time >= 0 & Time <= 15000),
    aes(x = Time, y = area, group = sample)
  ) +
    geom_line() +
    coord_cartesian(xlim = c(0, 15000)) +
    xlab("Time") +
    ylab("Area")

# Calculate growth rate taking into account initial tip variation
gr.res <- dplyr::mutate(gr.res, gr_per_tip = dplyr::if_else(tips.at.limit1 > 0, gr / tips.at.limit1, NA_real_))

# Plot growth rate taking into account initial tip variation
ggplot(gr.res, aes(x = tips.at.limit1, y = gr_per_tip)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    x = "Initial number of hyphal tips",
    y = "Growth rate per initial tip (area/hour/tip)",
    title = "Growth rate normalized by initial hyphal tip number"
  ) +
  theme_classic()

