# MetaGOF
A Bayesian Goodness-of-fit test for meta-analysis of rare events

## Access our method

### Introduction of files

***Main function***

MetaGOF_code: A main function to access our method, namely improved Pivotal quantities method (IPQ).

- We also provide three more frequentist-based approaches from Chen et al., 2015, C.-C. Wang and W.-C. Lee, 2019a and Naive method that employs Shaprio-Wilk test to the estimated log odds ratio directly.

***Stan files***:

- IW.stan: a stan file of using the *Inverse-Wishart* covariance prior.

- HIW.stan: a stan file of using the *Huang's Inverse-Wishart* covariance prior.

- SIW.stan: a stan file of using the *Scaled Inverse-Wishart* covariance prior.

- IND.stan: a stan file of using the *Indepedent* covariance prior.

***Data files***

- dat.bourassa1996.rda: Orignial meta-analysis of investigating the the hand-eye association conducted by Bourassa (1996), which contained 54 independent studies.
- eye_hand.data.rda: Cleaned data from dat.bourassa1996.rda. 
- t2d.data.rda: A meta-analysis of 20 independent studies to investigate the association between Type 2 diabetes mellitus and gestational diabetes conducted by Bellamy et al. (2009).

### How to use?

- Step 1: Download all files needed (at least stan files and the main function).
- Step 2: Install rstan, mvtnorm packages.
- Step 3: Open MetaGOF_code and run, where we recommend to use "IPQ" method with IND prior.

## References

- Johnson, V. E. (2007a). Bayesian model assessment using pivotal quantities. Bayesian Analysis, 2 (4).
- Yuan, Y., & Johnson, V. E. (2011). Goodness-of-fit diagnostics for bayesian hierarchical models. Biometrics, 68 (1), 156–164.
- Chen, Z., Zhang, G., & Li, J. (2015). Goodness-of-fit test for meta-analysis. Scientific Reports, 5 (1).
- Wang, C.-C., & Lee, W.-C. (2019a). Evaluation of the normality assumption in meta-analyses. American Journal of Epidemiology, 189 (3), 235–242.
- Bourassa, D. (1996). Handedness and eye-dominance: A meta-analysis of their relationship. Laterality, 1 (1), 5–34.
- Bellamy, L., Casas, J.-P., Hingorani, A. D., & Williams, D. (2009). Type 2 diabetes mellitus after gestational diabetes: A systematic review and meta-analysis. The Lancet, 373 (9677), 1773–1779.
- Liu, Y., & Xie, J. (2019). Cauchy combination test: A powerful test with analytic p-value calculation under arbitrary dependency structures. Journal of the American Statistical Association, 115 (529), 393–402.
