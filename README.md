# Conditional-Fixed-Income-Correlations-Applying-Extreme-Value-Theory

This repository contains the code of the "Conditional Fixed Income Correlations: Skews and Straddles" article that was published in the Journal of Fixed Income. The paper is accessible through this link:

https://jfi.pm-research.com/content/30/3/83

There are three R scripts that estimate conditional correlations using extreme value theory

- R_asymmetricCorrelations_negative_exceedances_Gumbel.R
- R_asymmetricCorrelations_positve_exceedances_Gumbel.R. 
- R_asymmetricCorrelations_cross_exceedances_Gumbel.R

We provide a truncated data file of monthly excess returns by the sectors and per unit of risk measured as spread duration. We only display three months of observations for copyright reasons. The data file is R_returnDataPerUnitRisk.csv. A more detailed discussion of the data can be found in the journal article.

# Conditional Extreme-Value Fixed Income Correlations Using Gumbel's Copula Specification

The three R scripts listed above estimate the maximum likelihood functions of excess returns beyond extremal thresholds by pre-specified pairs of Fixed Income assets in terms of their respective excess returns by unit of duration. The likelihood functions correspond to the Gumbel specification of the pairwise copulas. The extremal thresholds are parameterized to be the 65th, 75th, 85th and 90th percentiles. Each of the three R scripts are self-contained and can be run independently of each other. They only differ from one another in regards to the likelhood sub-function reflective of whether the analysis focuses on positive extremal co-movements in excess returns (R_asymmetricCorrelations_positve_exceedances_Gumbel.R. ), negative extremal co-movements in excess returns (R_asymmetricCorrelations_negative_exceedances_Gumbel.R) or cross-movements in extremal excess returns where one assets strongly outperforms and the other assets simultanously strongly underpreforms (R_asymmetricCorrelations_cross_exceedances_Gumbel.R). 

We also estimated alternative model specifications to estimate extremal co-movements of the Fixed Income assets under considerations. In particular, we estimated the Clayton and the Frank copula models. The scripts are available upon request.
