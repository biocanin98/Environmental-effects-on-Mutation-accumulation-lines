# Environmental-effects-on-Mutation-accumulation-lines

This project analyses environmental effects on Mutation accumulation lines.

For more information: https://pubmed.ncbi.nlm.nih.gov/36922001/

## Summary of the data, growth rate calculation and quality control check

1. For generating raw growth rate curves, use script oCelloscope_functions.R.

Note: different oCelloscope settings are used e.g. testdata <- load.ocelloscope.data2("Test1_scan_areas.csv", sep = ",", dec = ".", repetitions = 60, features = c("EstimatedVolume", "HyphaeArea", "NumberOfTips", "TipsPerArea", "TimestampInSeconds"))`

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

   Results in subfolder: Summary_Results


## Statistical analysis - Bayesian environmental-effect and invariance analysis for oCelloscope MA-line summaries

Place both MA_environmental_invariance_Bayesian.R and summary_results.zip in the
same folder in RStudio. Run source("MA_environmental_invariance_Bayesian.R").

Note: Requires brms, ggplot2, dplyr, tidyr, posterior, loo, and CmdStan/cmdstanr
or rstan installed in R

The script imports all seven environment folders from summary_results.zip;
it rechecks each well's adjusted R2 > 0.90, and combines only QC-passing
technical wells within line x environment x measurement x R1/R2 labels.
Uses the previously computed growth / initial tip with Student-t likelihood.

CRITICAL: anc_A has ZERO QC-passing wells at pH 8, and anc_a has only ONE.
No pH8 effects for mating-type A are estimable directly. The script
excludes these unsupported comparisons. pH8 mating-type a is very fragile.

## Results

1. Distribution of mutational effects: 12_environment_vs_control_posterior.csv file

   









