# Environmental-effects-on-Mutation-accumulation-lines

This project analyses environmental effects on Mutation accumulation lines. 

Research question: are the effects of parental-MA lines invariant in respect to their environment.

Stressful environments include: high temperature (38C and 40C), hypersalinity (0.5 M and 1 M of NaCl), pH (7 and 8)

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


## Statistical analysis - Bayesian environmental-effect and invariance analysis for oCelloscope MA-line summary

###Description

For this analysis I used Bayesian hierarchical models, using brms and Stan packages.

Growth rates were normalized by the initial number of hyphal tips. Only wells passing quality control (adjusted R² > 0.90) were retained, and technical replicates were averaged before fitting.

Two competing models were fitted using Student-t likelihood:

1. Invariant-effects model: assumes MA-lines effect remains constant across different environments, while accounting for environmental conditions, ancestral genetic background and measurement-level variation.

2. Environment-dependent model: allows MA-line effects to vary across environments, capturing genotype-by-environment interactions.

Predictive performance was evaluated using Pareto-smoothed importance sampling leave-one-out cross-validation (PSIS-LOO).



###Script
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

## Preliminary results

See: Results/Statistical_analysis_results/

1. Quality control check pass varies considerably between different environments (File: 01_QC_by_environment.csv Figure: 14_QC_pass_rates.png)

2. Mutational-effect distribution differs between environments, particularly at 38C (File: 11_environment_effect_distributions.csv, Figure: 15_environment_effect_distributions.png)

3. Significant increase of variance in mutational effects at temp 38C, while mean mutational effects show no systematic changes. (File: 12_environment_vs_control_posterior.csv, Figure: 15_environment_effect_distributions.png)

4. Four MA-lines with environment-dependent effect at 38C (L30, L35, L36 and L40) (File: 13_line_by_environment_invariance.csv, Figure: 17_line_effect_changes_95CI.png)

5. No universal environmental invariance - (File: 13_line_by_environment_invariance.csv, Figure: 16_invariance_probability_heatmap.png).

6. Mutational effects differ between ancestral genetic backgrounds at 38C (File: 10_line_effects_by_environment.csv, Figures: 15_environment_effect_distributions.png, 18_control_vs_environment_scatter.png)

7. The environment-dependent model showed no meaningful predictive advantage over the invariant-effects model (File: 09_Predictive_LOO_comparison.txt)


## Sources


Generation of MA lines:

https://pubmed.ncbi.nlm.nih.gov/36922001/

Coding source:

https://github.com/ikron/oCelloscope/tree/master

McElreath's Statistical Rethinking




















