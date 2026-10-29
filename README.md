# Environmental-effects-on-Mutation-accumulation-lines

This project analyses environmental effects on Mutation accumulation lines.

For more information: https://pubmed.ncbi.nlm.nih.gov/36922001/

## Data analysis

1. For generating raw growth rate curves, use script oCelloscope_functions.R.

For more information: https://github.com/ikron/oCelloscope/tree/master

2. Calculate growth rate taking into account variation of initial tip number - script: Growth_rate_calculations.R

Isolate windowed area

  `ggplot(
    dplyr::filter(longdata, Time >= 0 & Time <= 15000),
    aes(x = Time, y = area, group = sample)
  ) +
    geom_line() +
    coord_cartesian(xlim = c(0, 15000)) +
    xlab("Time") +
    ylab("Area")`

Calculate growth rate taking into account initial tip variation

`gr.res <- dplyr::mutate(gr.res, gr_per_tip = dplyr::if_else(tips.at.limit1 > 0, gr / tips.at.limit1, NA_real_))`

Plot growth rate taking into account initial tip variation

`ggplot(gr.res, aes(x = tips.at.limit1, y = gr_per_tip)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    x = "Initial number of hyphal tips",
    y = "Growth rate per initial tip (area/hour/tip)",
    title = "Growth rate normalized by initial hyphal tip number"
  ) +
  theme_classic()`

3. Quality control check (whether number of initial tips is sufficient) - script: Quality_control.R

4. Detailed comparison and quality control between two different measurements: compare_MA_measurements.R







