# Statistical Inference for the L1-Penalized Support Vector Machine in High Dimensions

## Abstract 

We develop a statistical inference framework for the $L_1$-penalized support vector machine (SVM) in the high-dimensional regime where both the sample size $n$ and the number of features $p$ diverge with a fixed ratio $\alpha = n/p$. Extending the general classification inference framework of Huang and Zeng (2025) to accommodate the non-differentiable hinge loss, we identify the required subgradient via the dual variables of a linear programming reformulation of the SVM. This characterization enables the construction of a debiased estimator whose components are asymptotically normal, facilitating valid confidence intervals and hypothesis tests for individual features. Extensive simulations across a range of covariance structures confirm the accuracy of the asymptotic theory. An application to a high-dimensional breast cancer gene expression dataset further demonstrates the method's ability to achieve simultaneous classification and variable selection in practice.

