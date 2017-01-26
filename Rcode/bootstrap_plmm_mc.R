# this script includes function to do bootstrapping for PLMMs
# it then simulated 3 datasets, and shows that the p values from bootstrapping are 
# higher than those from PLMMs (so to control type I errors).

library(geiger) # will load ape pkg automatically
library(phytools)
library(dplyr)
library(pez)
library(parallel)


#' The main purpose of this function is to do bootstrap for plmm
#' z.phy: a communityPGLMM object
#' vphy: the standardized var-cov phylogenetic matrix of species
#' n_simulation: how many bootstrap do you want?
#' NOTE: this function is time consuming.

bootstrap_z_plmm = function(z.phy, vphy, nsim = 1000, mc.core = 2, 
                            reltol = 1e-8, maxit = 1000){
  z_plmm = unname(z.phy$B.zscore["envi1:trait1",])
  
  # refit the data with a null model (i.e. no interaction term)
  z.no.interaction = communityPGLMM(abund ~ 1 + envi1 + trait1,
                                    data = z.phy$data, family = z.phy$family,
                                    sp = z.phy$data$sp, site = z.phy$data$site,
                                    random.effects = z.phy$random.effects,
                                    REML = z.phy$REML, verbose = F,
                                    reltol = reltol, maxit = maxit)
  
  nspp = nlevels(z.phy$sp)
  nsite = nlevels(z.phy$site)
  
  # SD of random terms: intercept, intercept phy, envi slopes, slope phy, site
  sd_re = sqrt(z.no.interaction$s2r)
  
  # simulate data sets
  sim_null_data = function(z.no.interaction = z.no.interaction){
    # the null model: 
    # Yi = alpha + aspp[i] + apspp[i]+ bsite[i] +
    # (beta1 + cspp[i] + cpspp[i]) env1site[i] + 
    # beta3 trait1spp[i] + ei 
    
    # simulate alpha + aspp[i] + apspp[i]+ bsite[i]
    a_spp = rnorm(n = nspp, mean = 0, sd = sd_re[1])
    a_spp_phy = MASS::mvrnorm(n = 1, mu = rep(0, nspp), Sigma = sd_re[2]^2 * vphy)
    mu_spp = z.no.interaction$B["(Intercept)", ] + a_spp + a_spp_phy # mean freq of sp
    
    b_site = rnorm(n = nsite, mean = 0, sd = sd_re[5])
    mu_i = rep(mu_spp, nsite) + rep(b_site, each = nspp) # each sp at each site
    # mu_i = rnorm(nspp * nsite, mean = mu_spp_site, sd = z.no.interaction$B.se[1]) # include SD of intercept
    
    # simulate (beta1 + cspp[i] + cpspp[i]) env1site[i]
    c_slope = rnorm(n = nspp, mean = 0, sd = sd_re[3])
    c_slope_phy = MASS::mvrnorm(n = 1, mu = rep(0, nspp), Sigma = sd_re[4]^4 * vphy)
    mu_envi_i = rep(z.no.interaction$B["envi1", ] + c_slope + c_slope_phy, nsite) * z.phy$data$envi1
    
    # simulate beta3 trait1spp[i]
    mu_trait_i = z.no.interaction$B["trait1", ] * z.phy$data$trait1
    
    # simulate ei
    e_resid = rnorm(nspp * nsite, mean = 0, sd = sqrt(z.no.interaction$s2resid))
    
    # put together
    y_i = mu_i + mu_envi_i + mu_trait_i + e_resid
    y_i
  }
  # end  sim_null_data
  
  z_bootstrap = mclapply(1:nsim, mc.cores = mc.core, function(x){
    set.seed(x)
    # simulate one data set
    y_sim = sim_null_data(z.no.interaction)
    dat = cbind(z.phy$data, y_sim)
    # then refit the data with full model
    fm = as.formula("y_sim ~ 1 + envi1 + trait1 + envi1 * trait1")
    z.interaction = communityPGLMM(fm,
                                   data = dat, family = z.phy$family,
                                   sp = z.phy$data$sp, site = z.phy$data$site,
                                   random.effects = z.phy$random.effects,
                                   REML = z.phy$REML, verbose = F,
                                   reltol = reltol, maxit = maxit, s2.init = z.phy$s2r)
    data.frame(z_refit = unname(z.interaction$B.zscore["envi1:trait1",]),
               p_refit = unname(z.interaction$B.pvalue["envi1:trait1",]))
  })
  z_bootstrap = do.call(rbind, args = z_bootstrap)
  
  n_larger = sum(abs(z_bootstrap$z_refit) > abs(z_plmm))
  p = n_larger/(nsim + 1)
  # return these statistics
  list(statistics = data.frame(z_plmm = z_plmm, 
             p_plmm = unname(z.phy$B.pvalue["envi1:trait1",]),
             p_bootstrap = p,
             mean_z_null = mean(z_bootstrap$z_refit, na.rm = T),
             sd_z_null = sd(z_bootstrap$z_refit, na.rm = T),
             mean_p_null = mean(z_bootstrap$p_refit, na.rm = T)),
       z_bootstrap = z_bootstrap)
}

# examples ----

# This function to simulate a dataset, and fit a pglmm to it
sim_example = function(nspp = 10, nsite = 10, physig.trait1 = TRUE,
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
  
  z.phy = communityPGLMM(abund ~ 1 + envi1 + trait1 + trait1 * envi1,
                         data = dat, family = family,
                         sp = dat$sp, site = dat$site,
                         random.effects = list(re.sp, re.sp.phy, re.envi1,
                                               re.envi1.phy, re.site),
                         REML = REML, verbose = F,
                         reltol = reltol, maxit = maxit)
  
  list(z.phy, vphy)
}

# do not run it if you cannot wait...
# n.core = 20 # change it for your purpose
# z.phy = sim_example(nspp = 10, nsite = 10)
# sp10site10_n1000 = bootstrap_z_plmm(z.phy[[1]], z.phy[[2]], nsim = 1000, mc.core = n.core)
# if(!dir.exists("../rds_bootstrap")) dir.create("../rds_bootstrap")
# saveRDS(sp10site10_n1000, file = "../rds_bootstrap/sp10site10_n1000.rds")
# z.phy = sim_example(nspp = 20, nsite = 20)
# sp10site10_n1000 = bootstrap_z_plmm(z.phy[[1]], z.phy[[2]], nsim = 1000, mc.core = n.core)
# saveRDS(sp10site10_n1000, file = "../rds_bootstrap/sp20site20_n1000.rds")
# z.phy = sim_example(nspp = 30, nsite = 30)
# sp10site10_n1000 = bootstrap_z_plmm(z.phy[[1]], z.phy[[2]], nsim = 1000, mc.core = n.core)
# saveRDS(sp10site10_n1000, file = "../rds_bootstrap/sp30site30_n1000.rds")

sp10site10_n1000 = readRDS(file = "../rds_bootstrap/sp10site10_n1000.rds")
sp10site10_n1000[[1]]
#      z_plmm    p_plmm p_bootstrap   mean_z_null sd_z_null mean_p_null
#  0.02391723 0.9809186    0.984016 -0.01404447  1.395022    0.431283
mean(sp10site10_n1000$z_bootstrap$p_refit < 0.05)
# 0.136
hist(sp10site10_n1000$z_bootstrap$p_refit)
x = seq(-10, 10, length.out = 500)
y = dnorm(x, mean = 0, sd = 1)
hist(sp10site10_n1000$z_bootstrap$z_refit, freq = F, ylim = c(0, 1))
lines(x, y, col = "blue")

sp20site20_n1000 = readRDS(file = "../rds_bootstrap/sp20site20_n1000.rds")
sp20site20_n1000[[1]]
#      z_plmm    p_plmm p_bootstrap   mean_z_null sd_z_null mean_p_null
#  -0.3125529 0.7546204   0.8041958   0.0355099   1.25491   0.4412774
mean(sp20site20_n1000$z_bootstrap$p_refit < 0.05)
# 0.115
hist(sp20site20_n1000$z_bootstrap$p_refit)


sp30site30_n1000 = readRDS(file = "../rds_bootstrap/sp30site30_n1000.rds")
sp30site30_n1000[[1]]
#      z_plmm    p_plmm  p_bootstrap   mean_z_null sd_z_null mean_p_null
# 1 0.2234509 0.8231846   0.8421578 -0.02165727  1.140824   0.4670882
mean(sp30site30_n1000$z_bootstrap$p_refit < 0.05)
# 0.084
hist(sp30site30_n1000$z_bootstrap$p_refit)

par(mfrow = c(2,2))
hist(sp10site10_n1000$z_bootstrap$z_refit, freq = F, ylim = c(0, 0.5), breaks = 100)
lines(x, y, col = "blue")
hist(sp20site20_n1000$z_bootstrap$z_refit, freq = F, ylim = c(0, 0.5), breaks = 100)
lines(x, y, col = "blue")
hist(sp30site30_n1000$z_bootstrap$z_refit, freq = F, ylim = c(0, 0.5), breaks = 100)
lines(x, y, col = "blue")
# dev.copy2pdf(file = "bootstrap_z.pdf")
par(mfrow = c(1,1))

rbind(sp10site10_n1000[[1]], 
      sp20site20_n1000[[1]],
      sp30site30_n1000[[1]]) %>% 
  mutate(sim = c("sp10site10", "sp20site20", "sp30site30"))
#        z_plmm    p_plmm p_bootstrap mean_z_null sd_z_null mean_p_null        sim
# 1  0.02391723 0.9809186   0.9840160 -0.01404447  1.395022   0.4312830 sp10site10
# 2 -0.31255288 0.7546204   0.8041958  0.03550990  1.254910   0.4412774 sp20site20
# 3  0.22345086 0.8231846   0.8421578 -0.02165727  1.140824   0.4670882 sp30site30
