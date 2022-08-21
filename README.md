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

### How to use our method?

- Step 1: Download all files needed (at least stan files and the main function).
- Step 2: Install rstan, mvtnorm packages.
- Step 3: Open MetaGOF_code and run, where we recommend to use "IPQ" method with IND prior.

