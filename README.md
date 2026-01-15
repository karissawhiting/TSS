# Two-stage Subsampling Variable Selection 

This code implements the two-stage high-dimensional variable selection method (TSS) for generalized linear models. The user needs to supply the high-dimensional matrix of predictors x, the outcome y, the number of subsamples used in each stage, the cutoff used for the selection threshold in Stage 1 and the distribution of the outcome. 


# Notes on Needed Code Updates

- Need to Check for 0 variance before `plsRglm`
- Need to check that scad returned results before `plsRglm`
- Update threshold logic (sometimes we use top 10)



Capanu, M., Giurcanu, M., Begg, C. B., & Gönen, M. (2025). Two-stage subsampling variable selection for sparse high-dimensional generalized linear models. *Statistical Methods in Medical Research, 34*(7), 1504–1521. https://doi.org/10.1177/09622802251343597
 