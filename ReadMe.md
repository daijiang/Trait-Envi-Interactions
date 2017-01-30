This repository includes R scripts to simulate datasets, fit models (LMMs and PLMMs), and plot the results. It also includes R functions to do bootstrapping on PLMMs and examples about how to use it.

Theses scripts are under MIT license and come with no guaranntees. Please cite the paper if you use the code elsewhere.

All questions should be directed to Daijiang Li (daijianglee@gmail.com).

- `Rcode` folder: all R codes used in the main text.
    + `simulation.R`: simulates datasets, fits models (LMMs and PLMMs), and saves results.
    + `plot.R`: plot the results, reproduce all figures from the paper.
    + `bootstrap_plmm_mc.R`: do bootstrapping on PLMMs, with multiple cores as an option.
- `rds` folder: model fitting results from LMMs and PLMMs with varing number of species (i.e. Type I errors).
- `rds2` folder: model fitting results from LMMs and PLMMs with varing true parameters (i.e. Type II errors).
- `rds_bootstrap` folder: results of bootstrapping on PLMMs.
- `appendixS1` folder: R code and data to show that they way to simulate the phylogeny does not affect our results and conclusions. The R codes are mostly the same as those in the main text, with their ways to simulate phylogeny changed.

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


MIT License

Copyright (c) [2017] [Daijiang Li]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
