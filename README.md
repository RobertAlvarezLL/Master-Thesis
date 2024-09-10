# Master’s Thesis: Time-Varying Effects of Extreme Weather Shocks on Main Macroeconomic Indicators of Spain

## Code Overview

The code provided in this repository illustrates our approach to implementing Structural Vector Autoregression (SVAR) models with sign restrictions in R, addressing the gap left by the absence of a pre-existing R package that fully accommodates our analytical needs. This custom approach combines the use of the `vars` package for initial VAR estimation with custom functions to impose sign restrictions and identify valid structural matrices. The iterative process ensures that the imposed economic theory-based constraints are met.

For the Time-Varying Parameters (TVP) Bayesian VAR model, the code consists of the implementation in the **BVARS** package (Krueger, 2015) of the Primiceri (2005) model, including the corrigendum by Del Negro and Primiceri (2015). Additionally, a modification to the Primiceri model has been implemented following the method proposed by Arias, Rubio-Ramírez, and Waggoner (2014) and further applied by Leiva-Leon and Uzeda (2023) to incorporate sign restrictions into the TVP-VAR framework. This adaptation allows us to analyze the evolving impact of extreme weather shocks on macroeconomic variables in a dynamic and flexible manner.




