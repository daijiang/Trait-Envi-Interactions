# This script will simulate datasets and then fit LMMs and PLMMs.
# Warnning: this script is really time-consuming if you want to re-run it.
# As we setted seeds in simulations, you can reproduce all results. 
# But it may take you weeks to finish the simulations.
# We have saved our model fitting results in the rds and rds2 folders if you are interested.

library(geiger) # will load ape pkg automatically
library(phytools)
library(plyr)
library(dplyr)
# library(ggplot2)
library(tidyr)
library(pez)
library(parallel)

## this function will get an initial set of s2 for random effects
## for following simulations, to improve speed.
sim_init = function(nspp = 30, nsite = 30, physig.trait1 = TRUE,
                    physig.trait2 = TRUE,
                    alpha = 1, beta_1 = 1, beta_2 = 1, beta_3 = 1,
                    beta_4 = 1, beta_5 = 0, beta_6 = 1,
                    family = "gaussian", reltol = 10^-8,
                    maxit = 1000, REML = FALSE){
  # simulate the phylogeny
  phy = sim.bdtree(b = 1, d = 0, stop = "taxa", n = nspp, extinct = FALSE)
  phy$tip.label = paste0("sp", 1:nspp)
  vphy = vcv(phy) # var-cov matrix
  vphy = vphy/(det(vphy)^(1/nspp)) # to standardize the matrix
  corphy = cov2cor(vphy) # corr matrix
  # simulate traits
  if(physig.trait1){
    trait1 = fastBM(tree = phy, mu = 0, sig2 = 1)
  } else {
    trait1 = rnorm(nspp, 0, 1)
  }
  if(physig.trait2){
    trait2 = fastBM(tree = phy, mu = 0, sig2 = 1)
  } else {
    trait2 = rnorm(nspp, 0, 1)
  }
  trait = data.frame(sp = paste0("sp", 1:nspp),
                     trait1 = trait1,
                     trait2 = trait2, stringsAsFactors = FALSE)
  
  # simulate two envi variables
  envi = data.frame(site = paste0("site", 1:nsite),
                    envi1 = sort(runif(n = nsite, min = -1, max = 1)),
                    envi2 = rnorm(n = nsite),
                    stringsAsFactors = FALSE)
  
  # simulate abundance data
  dat = expand.grid(sp = paste0("sp", 1:nspp),
                    site = paste0("site", 1:nsite)) %>%
    mutate(sp = as.character(sp), site = as.character(site))
  dat = left_join(dat, trait, by = "sp") %>%
    left_join(envi, by = "site")
  dat = mutate(dat,
               abund = alpha + beta_1 * envi1 + beta_2 * envi2 +
                 beta_3 * trait1  + beta_4 * trait2 +
                 beta_5 * envi1 * trait1 +
                 beta_6 * envi1 * trait2 +
                 rnorm(nspp * nsite, 0, 1)) %>%
    mutate(sp = factor(sp, levels = paste0("sp", 1:nspp)),
           site = factor(site, levels = paste0("site", 1:nsite)))
  
  # model fitting
  re.site <- list(1, site = dat$site, covar = diag(nsite))
  re.sp <- list(1, sp = dat$sp, covar = diag(nspp))
  re.sp.phy <- list(1, sp = dat$sp, covar = vphy)
  re.envi1 = list(dat$envi1, sp = dat$sp, covar = diag(nspp))
  re.envi1.phy <- list(dat$envi1, sp = dat$sp, covar = vphy)
  z.no.phy = communityPGLMM(abund ~ 1 + envi1 + trait1 + trait1 * envi1,
                            data = dat, family = family,
                            sp = dat$sp, site = dat$site,
                            random.effects = list(re.sp, re.envi1, re.site),
                            REML = REML, verbose = F,
                            reltol = reltol, maxit = maxit)[c("B", "B.se",
                                                              "B.zscore",
                                                              "B.pvalue",
                                                              "s2r", "s2resid")]
  z.phy = communityPGLMM(abund ~ 1 + envi1 + trait1 + trait1 * envi1,
                         data = dat, family = family,
                         sp = dat$sp, site = dat$site,
                         random.effects = list(re.sp, re.sp.phy, re.envi1,
                                               re.envi1.phy, re.site),
                         REML = REML, verbose = F,
                         reltol = reltol, maxit = maxit)[c("B", "B.se",
                                                           "B.zscore",
                                                           "B.pvalue",
                                                           "s2r", "s2resid")]
  
  list(z.no.phy$s2r, z.phy$s2r)
}

## This function does one simulation,
## which will be nested with the sim_n() function.
sim_1 = function(nspp = 30, nsite = 30, physig.trait1 = TRUE,
                 physig.trait2 = TRUE,
                 alpha = 1, beta_1 = 1, beta_2 = 1, beta_3 = 1,
                 beta_4 = 1, beta_5 = 0, beta_6 = 1,
                 family = "gaussian", reltol = 10^-8,
                 maxit = 1000, REML = FALSE, s2_init_no_phy, s2_init_phy){
  # simulate the phylogeny
  phy = sim.bdtree(b = 1, d = 0, stop = "taxa", n = nspp, extinct = FALSE)
  phy$tip.label = paste0("sp", 1:nspp)
  vphy = vcv(phy) # var-cov matrix
  vphy = vphy/(det(vphy)^(1/nspp)) # to standardize the matrix
  corphy = cov2cor(vphy) # corr matrix
  # simulate traits
  if(physig.trait1){
    trait1 = fastBM(tree = phy, mu = 0, sig2 = 1)
  } else {
    trait1 = rnorm(nspp, 0, 1)
  }
  if(physig.trait2){
    trait2 = fastBM(tree = phy, mu = 0, sig2 = 1)
  } else {
    trait2 = rnorm(nspp, 0, 1)
  }
  trait = data.frame(sp = paste0("sp", 1:nspp),
                     trait1 = trait1,
                     trait2 = trait2, stringsAsFactors = FALSE)
  
  # simulate two envi variables
  envi = data.frame(site = paste0("site", 1:nsite),
                    envi1 = sort(runif(n = nsite, min = -1, max = 1)),
                    envi2 = rnorm(n = nsite),
                    stringsAsFactors = FALSE)
  
  # simulate abundance data
  dat = expand.grid(sp = paste0("sp", 1:nspp), site = paste0("site", 1:nsite)) %>%
    mutate(sp = as.character(sp), site = as.character(site))
  dat = left_join(dat, trait, by = "sp") %>%
    left_join(envi, by = "site")
  dat = mutate(dat,
               abund = alpha + beta_1 * envi1 + beta_2 * envi2 +
                 beta_3 * trait1  + beta_4 * trait2 +
                 beta_5 * envi1 * trait1 +
                 beta_6 * envi1 * trait2 +
                 rnorm(nspp * nsite, 0, 1)) %>%
    mutate(sp = factor(sp, levels = paste0("sp", 1:nspp)),
           site = factor(site, levels = paste0("site", 1:nsite)))
  
  # model fitting
  re.site <- list(1, site = dat$site, covar = diag(nsite))
  re.sp <- list(1, sp = dat$sp, covar = diag(nspp))
  re.sp.phy <- list(1, sp = dat$sp, covar = vphy)
  re.envi1 = list(dat$envi1, sp = dat$sp, covar = diag(nspp))
  re.envi1.phy <- list(dat$envi1, sp = dat$sp, covar = vphy)
  z.no.phy = communityPGLMM(abund ~ 1 + envi1 + trait1 + trait1 * envi1,
                            data = dat, family = family,
                            sp = dat$sp, site = dat$site,
                            random.effects = list(re.sp, re.envi1, re.site),
                            REML = REML, verbose = F,
                            reltol = reltol, maxit = maxit,
                            s2.init = s2_init_no_phy)[c("B", "B.se",
                                                        "B.zscore",
                                                        "B.pvalue",
                                                        "s2r", "s2resid")]
  z.phy = communityPGLMM(abund ~ 1 + envi1 + trait1 + trait1 * envi1,
                         data = dat, family = family,
                         sp = dat$sp, site = dat$site,
                         random.effects = list(re.sp, re.sp.phy, re.envi1,
                                               re.envi1.phy, re.site),
                         REML = REML, verbose = F,
                         reltol = reltol, maxit = maxit,
                         s2.init = s2_init_phy)[c("B", "B.se",
                                                  "B.zscore",
                                                  "B.pvalue",
                                                  "s2r", "s2resid")]
  
  extract_coef_pglmm = function(pglmm){
    coef.df = data.frame(estimate = pglmm$B, SE = pglmm$B.se,
                         Z = pglmm$B.zscore, Pvalue = pglmm$B.pvalue)
    rownames(coef.df)[1] = "Intercept"
    coef.df$terms = rownames(coef.df)
    rownames(coef.df) = NULL
    coef.df[, c(5,1:4)] # reorder columns
  }
  
#  extract random terms sd
  est.RE = rbind(data.frame(terms = c("a.RE", "env.RE", "site.RE", "resid.RE"),
             var = c(z.no.phy$s2r, z.no.phy$s2resid),
             trt = "no_phylo"),
    data.frame(terms = c("a.RE", "a.phy.RE", "env.RE", "env.phy.RE", "site.RE", "resid.RE"),
               var = c(z.phy$s2r, z.phy$s2resid),
               trt = "phylo"))

  
  est.FE = bind_rows(mutate(extract_coef_pglmm(z.no.phy), trt = "no_phylo"),
                     mutate(extract_coef_pglmm(z.phy), trt = "phylo")) %>%
    left_join(data.frame(terms = c("Intercept", "envi1", "trait1", "envi1:trait1"),
                         real.val = c(alpha, beta_1, beta_3, beta_5),
                         stringsAsFactors = F), by = "terms")

  list(est.fixed.eff = est.FE, nspp = nspp, nsite = nsite, est.random.eff = est.RE)
}

sim_n = function(n = 10, mc.cores = 4, nspp = 30, nsite = 30,
                 physig.trait1 = TRUE, physig.trait2 = TRUE,
                 alpha = 1, beta_1 = 1, beta_2 = 1, beta_3 = 1,
                 beta_4 = 1, beta_5 = 0, beta_6 = 1,
                 family = "gaussian", reltol = 10^-8,
                 maxit = 1000, REML = FALSE){
  # to get an initial value that will be used later
  s2_init = sim_init(nspp = nspp, nsite = nsite,
                     physig.trait1 = physig.trait1,
                     physig.trait2 = physig.trait2,
                     alpha = alpha, beta_1 = beta_1, beta_2 = beta_2,
                     beta_3 = beta_3, beta_4 = beta_4, beta_5 = beta_5,
                     beta_6 = beta_6, family = family, reltol = reltol,
                     maxit = maxit, REML = REML)
  # do multiple simulations (does not work on Windows OS)
  xt = mclapply(seq_len(n), function(x) {
    set.seed(x)
    sim_1(nspp = nspp, nsite = nsite,
          physig.trait1 = physig.trait1, physig.trait2 = physig.trait2,
          alpha = alpha, beta_1 = beta_1, beta_2 = beta_2,
          beta_3 = beta_3, beta_4 = beta_4, beta_5 = beta_5,
          beta_6 = beta_6, family = family, reltol = reltol,
          maxit = maxit, REML = REML, s2_init_no_phy = s2_init[[1]],
          s2_init_phy = s2_init[[2]])
  }, mc.cores = mc.cores)
  # to extract fixed effects
  est.FE = mutate(ldply(xt, function(x) x$est.fixed.eff),
                  item = rep(1:length(xt), each = 8))
  est.RE = mutate(ldply(xt, function(x) x$est.random.eff),
                  item = rep(1:length(xt), each = 10))
  output = list(est.fixed.eff = est.FE,
                nspp = xt[[1]]$nspp, 
                nsite = xt[[1]]$nsite, 
                est.random.eff = est.RE)
  
  # save output on disk
  # if beta_5 = 0, then this is for Type I errors
  if(beta_5 == 0){
    if(!dir.exists("./rds")) dir.create("rds")
    saveRDS(output, file = paste("rds/sp", nspp, "site", nsite, "_",
                                 as.integer(physig.trait1),
                                 as.integer(physig.trait2),
                                 ".rds", sep = ""))
  } else{ # this is for power analyses
    if(!dir.exists("./rds2")) dir.create("rds2")
    saveRDS(output, file = paste("rds2/sp", nspp, "site", nsite, "_",
                                 "beta_5", beta_5, "_", as.integer(physig.trait1),
                                 as.integer(physig.trait2),
                                 ".rds", sep = ""))
  }
  output
}

n_cpus = 18
n_sim = 1000
# type I error
sim_n(n = n_sim, nspp = 20, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 20, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 20, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 20, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = F)

sim_n(n = n_sim, nspp = 30, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 30, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 30, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 30, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = F)

sim_n(n = n_sim, nspp = 40, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 40, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 40, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 40, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = F)

sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = F)
# these four simulations will be used for type II error above.

sim_n(n = n_sim, nspp = 60, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 60, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 60, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 60, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = F)

sim_n(n = n_sim, nspp = 70, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 70, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 70, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 70, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = F)

sim_n(n = n_sim, nspp = 80, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 80, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 80, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 80, nsite = 30, mc.cores = n_cpus, beta_5 = 0, physig.trait1 = F, physig.trait2 = F)

# type II error
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.25, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.25, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.25, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.25, physig.trait1 = F, physig.trait2 = F)

sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.50, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.50, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.50, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.50, physig.trait1 = F, physig.trait2 = F)

sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.75, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.75, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.75, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 0.75, physig.trait1 = F, physig.trait2 = F)

sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 1.00, physig.trait1 = T, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 1.00, physig.trait1 = F, physig.trait2 = T)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 1.00, physig.trait1 = T, physig.trait2 = F)
sim_n(n = n_sim, nspp = 50, nsite = 30, mc.cores = n_cpus, beta_5 = 1.00, physig.trait1 = F, physig.trait2 = F)