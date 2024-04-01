PPMx model write up 

# Context
## Complex traits
Complex traits have causal linkages dispesed across many locations on the genome (as opposed to a few or a single). Height is a commonly cited complex trait, whereas many of the first studied monogenic traits relate to functionality of hemoglobin.

## Approaches in predicting complex traits from genotyping data
Genotyping data regularly contain sequencing information on 10k - 1M+ single nucleotide polymorphisms (SNPs). These studies are generally made up of <10k subjects leading to problems in estimability. To get around this methods can either "prune" or utilize the correlated structure of the SNPs, called Linkage disequilibrium (LD) to create more reliable phenootype estimates called Polygenic Risk Scores (PRS). They generally follow theses steps:

1. Estimate marginal effects of SNPs (can be from different study or same study with sample splitting)
$$
  ∀ j \in \{1, \dots, m\} \hat β^{GWAS} argmin_{β_j} (Y - X^{train}_{.j} β_j)^2
$$

2. (Prune) SNP selection via effect size or p-value thresholds

3. Construct a PRS (Can be more complicated when accounting for LD) 
$$
  PRS = X \hat \beta^{GWAS}
$$

### Difficulties with heterogeneity
Heterogeneous populations can lead to difficulties in constructing a PRS and can arise from the following conditions 
- Different allele frequencies $(X)$ as occurs in different ancestries
  - Adjust for PC's of X
  - Standardize X
- Different SNP effects 
  - Reestimating the SNP effects assuming that SNP effects will be similar (VIPRS)

Both of these can change the association between a PRS and a phenotype.


## Proposed approach (PPMx)
Definitions
- Bayesian - Statistical method combining observed data with prior knowledge
- Non-parametric - The model contains too many parameters to perform inference making the model effectively non-parametric

Product Partition Model with covariates (PPMx) is a Bayesian nonparametric approach that combines clustering with regression. PPMx detects clusters based on the association between $Y$ and $X$.

Using this method is proposed to help with PRS modeling by
- Detecting differences in scaling ($X$) associated with 
- Detecting differences in association between genotype and phenotype ($β$)


## The PPMx model
- $ρ_n$ partition of clusters
- $\{S_k: k ∈ 1,...,K\}$ clusters of subjects
- $c(S_k)$ coherence function for cluster k. (often polya urn)
- $g(x_j^*)$ similarity function for covariates in cluster j

### Prior
Extends the polya urn prior
$$
  p(ρ_n = \{S_1, …, S_K\} | X) ∝ ∏_{j=1}^K c(S_j) g(x_j^*) \tag{prior} \\
  = M^K∏(|S_j|-1)! g(x_j^*)
$$

Notice that the first portion is the Polya Urn again. So we focus on the similarity function $g(x_j^*)$.

### Similarity function
Let $\tilde ξ$ be latent variable purely used for simplifcation by allowing us to model $x_i$ based on cluster specific parameter $ξ$ using $q(x|ξ)$.

$$
  g(x_j^*) = ∫∏_{i∈ S_j}q(x_i|\tilde ξ_j^*)q(\tilde ξ_j^*)d\tilde ξ_j^*
$$

Using conjugacy we can simplify this to

$$
  g(x_j^*) = \frac{ ∏_{i∈ S_j}q(x_i|\tilde ξ_j^*)q(\tilde ξ_j^*)}{q(\tilde ξ_j^*|x^*_j)}
$$


### Choice of $ξ$

$$
  ξ_j^* = (m_j, v_j) \\
  q(x_i | ξ_j^*) = N(x_i | m_j, v_j) \\
  q(m_j, v_j) = NInvWishart()
$$

### Full model statement
$$
  y_i | \alpha^*, \sigma^2 \sim N(\alpha_{S_i}x_i, \sigma^2) \\
  \alpha^*_j | \sigma^2 \sim N(a, \tau_0^2\sigma^2)\\
  \sigma^2 \sim IGamma(v_0, \lambda_0) \\
  p(\rho_n = \{S_1,...,S_K\} | X) \propto \prod_{j=1}^K c(S_j)g(x_j^*)
$$

## Problem statement and method extension
PPMx only incorporates variables that are informative of cluster formation. We would like to account for scanner effects in away that doesn't inform our clusters to make sure the clusters relate more biologically relevant information. So we propose to extend PPMx by adding another set of covariates which are non-cluster informative. If promising on a single phenotype, we would like to use this model to attempt mapping between imaging IDP's and clinical/ behavioral outcomes.

$$
  y_i | α^*, β_S, σ^2 \sim N(α_{S_i}x_i + β_sx_{S,i}, σ^2) \\
  β_S | η \sim N(0, η) \\
  η | c,d \sim IG(c,d) \\ 
  \alpha^*_j | \sigma^2 \sim N(a, \tau_0^2\sigma^2)\\
  \sigma^2 \sim IGamma(v_0, \lambda_0) \\
  p(\rho_n = \{S_1,...,S_K\} | X) \propto \prod_{j=1}^K c(S_j)g(x_j^*)
$$

Fitting process outlline, residualize out site/scanner effects first then fit PPMx model.

### Pros and Cons

------------------
| Pros | Cons |
| ---- | ---- |
| Computationally efficient. Allows for an MCMC solution  | Doesn't reestimate $β$|
| Fully Bayesian | Restrictive assumption. Assumes differing $β$ will result in a different slope between clusters |
| Detect meaningful clustering in a data driven way |  |


