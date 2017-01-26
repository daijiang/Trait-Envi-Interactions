This repository includes R scripts to simulate datasets, fit models (LMMs and PLMMs), and plot the results. It also includes R functions to do bootstrapping on PLMMs and examples about how to use it.

The code comes with no guaranntees and all questions should be directed to Daijiang Li (daijianglee@gmail.com).

Please cite the paper if you use the code elsewhere.

- `Rcode` folder: all R codes used.
    + `simulation.R`: simulates datasets, fits models (LMMs and PLMMs), and saves results.
    + `plot.R`: plot the results, reproduce all figures from the paper.
    + `bootstrap_plmm_mc.R`: do bootstrapping on PLMMs, with multiple cores as an option.
- `rds` folder: model fitting results from LMMs and PLMMs with varing number of species (i.e. Type I errors).
- `rds2` folder: model fitting results from LMMs and PLMMs with varing true parameters (i.e. Type II errors).
- `rds_bootstrap` folder: results of bootstrapping on PLMMs.

        `sessionInfo()`
        R version 3.3.2 (2016-10-31)
        Platform: x86_64-apple-darwin13.4.0 (64-bit)
        Running under: macOS Sierra 10.12.2
        
        locale:
        [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
        
        attached base packages:
        [1] stats     graphics  grDevices utils     datasets  methods   base     
        
        other attached packages:
        [1] viridis_0.3.4 cowplot_0.6.2 ggplot2_2.2.0 tidyr_0.6.0   dplyr_0.5.0   plyr_1.8.4   
        
        loaded via a namespace (and not attached):
         [1] Rcpp_0.12.8      digest_0.6.10    assertthat_0.1   MASS_7.3-45      grid_3.3.2       R6_2.2.0         gtable_0.2.0     DBI_0.5          magrittr_1.5    
        [10] scales_0.4.1     stringi_1.1.1    reshape2_1.4.1   lazyeval_0.2.0   labeling_0.3     tools_3.3.2      stringr_1.1.0    munsell_0.4.3    colorspace_1.2-6
        [19] gridExtra_2.2.1  tibble_1.2