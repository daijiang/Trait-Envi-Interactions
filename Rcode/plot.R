# this script will reproduce all figures in the paper.
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)

## read rds files from type I errors: varied nspp ----
rdss1 = list.files("../rds", full.names = T)
sim.FE1 = vector("list", length = length(rdss1))
for (i in 1:length(rdss1)){
  sim.FE1[[i]] = readRDS(file = rdss1[i])$est.fixed.eff %>%
    mutate(nspp = as.numeric(gsub("^.*sp([0-9]{2})site.*$", "\\1", rdss1[i])),
           nsite = as.numeric(gsub("^.*site([0-9]{2}).*$", "\\1", rdss1[i])),
           catg = gsub("^.*_([0-9]{2}).*$", "\\1", rdss1[i]))
}
sim.FE1 = tbl_df(ldply(sim.FE1)) %>% filter(terms == "envi1:trait1")

prop_better1 = sim.FE1 %>%
  select(estimate, trt, real.val, item, nspp, catg) %>%
  mutate(diff = abs(estimate - real.val)) %>%
  select(-estimate) %>%
  spread(key = trt, value = diff)
prop_better1 = adply(prop_better1, 1, function(x){# take a while
  if(x[5] > x[6]) better_trt = "phylo"
  if(x[5] < x[6]) better_trt = "no_phylo"
  if(x[5] == x[6]) better_trt = sample(c("phylo", "no_phylo"), 1)
  better_trt
})
prop_better1 = prop_better1 %>% tbl_df() %>%
  select(nspp, catg, V1) %>%
  group_by(nspp, catg) %>%
  summarise(PLMM = sum(V1 == "phylo")/max(sim.FE1$item),
            LMM = sum(V1 == "no_phylo")/max(sim.FE1$item)) %>%
  gather(key = "Models", value = "prop", PLMM, LMM) %>%
  ungroup() %>%
  mutate(nspp = as.character(nspp))
prop_better1$catg[prop_better1$catg == "00"] = "trait1: I; trait2: I"
prop_better1$catg[prop_better1$catg == "01"] = "trait1: I; trait2: C"
prop_better1$catg[prop_better1$catg == "10"] = "trait1: C; trait2: I"
prop_better1$catg[prop_better1$catg == "11"] = "trait1: C; trait2: C"
prop_better1$catg = factor(prop_better1$catg,
                           levels = c("trait1: I; trait2: I",
                                      "trait1: C; trait2: I",
                                      "trait1: I; trait2: C",
                                      "trait1: C; trait2: C"))
prop_better1$Models[prop_better1$Models == "PLMM"] = "PLMM better"
prop_better1$Models[prop_better1$Models == "LMM"] = "LMM better"
# data to plot type I error
type_errs = sim.FE1 %>% 
  group_by(trt, nspp, nsite, catg) %>%
  summarise(nitem = max(item),
            mid.est = median(estimate),
            ave.est = mean(estimate),
            ave.est.low = mean(estimate) - 1 * sd(estimate),
            ave.est.high = mean(estimate) + 1 * sd(estimate),
            reject.prop = sum(Pvalue <= 0.05)/nitem) %>% ungroup()
type_errs$catg2 = type_errs$catg
type_errs$catg[type_errs$catg == "00"] = "trait1: I; trait2: I"
type_errs$catg[type_errs$catg == "01"] = "trait1: I; trait2: C"
type_errs$catg[type_errs$catg == "10"] = "trait1: C; trait2: I"
type_errs$catg[type_errs$catg == "11"] = "trait1: C; trait2: C"
type_errs$catg = factor(type_errs$catg,
                        levels = c("trait1: I; trait2: I",
                                   "trait1: C; trait2: I",
                                   "trait1: I; trait2: C",
                                   "trait1: C; trait2: C"))
type_errs = mutate(type_errs, models = trt)
type_errs$models[type_errs$models == "no_phylo"] = "LMM"
type_errs$models[type_errs$models == "phylo"] = "PLMM"


## read rds files from power / type II errors: varied beta 5 ----
rdss2 = list.files("../rds2", full.names = T)
sim.FE2 = vector("list", length = length(rdss2))
for (i in 1:length(rdss2)){
  sim.FE2[[i]] = readRDS(file = rdss2[i])$est.fixed.eff %>%
    mutate(nspp = as.numeric(gsub("^.*sp([0-9]{2})site.*$", "\\1", rdss2[i])),
           nsite = as.numeric(gsub("^.*site([0-9]{2}).*$", "\\1", rdss2[i])),
           beta5 = as.numeric(gsub("^.*beta_5([0-9.]{1,4}).*$", "\\1", rdss2[i])),
           catg = gsub("^.*_([0-9]{2}).*$", "\\1", rdss2[i]))
}
sim.FE2 = tbl_df(ldply(sim.FE2)) %>% filter(terms == "envi1:trait1")

prop_better2 = sim.FE2 %>%
  select(estimate, trt, real.val, item, nspp, catg) %>%
  mutate(diff = abs(estimate - real.val)) %>%
  select(-estimate) %>%
  spread(key = trt, value = diff)
prop_better2 = adply(prop_better2, 1, function(x){
  if(x[5] > x[6]) better_trt = "phylo"
  if(x[5] < x[6]) better_trt = "no_phylo"
  if(x[5] == x[6]) better_trt = sample(c("phylo", "no_phylo"), 1)
  better_trt
})
prop_better2 = prop_better2 %>% tbl_df() %>%
  select(real.val, catg, V1) %>%
  group_by(real.val, catg) %>%
  summarise(PLMM = sum(V1 == "phylo")/max(sim.FE2$item),
            LMM = sum(V1 == "no_phylo")/max(sim.FE2$item)) %>%
  gather(key = "Models", value = "prop", PLMM, LMM) %>%
  ungroup() %>%
  mutate(real.val = as.character(real.val))
prop_better2$catg[prop_better2$catg == "00"] = "trait1: I; trait2: I"
prop_better2$catg[prop_better2$catg == "01"] = "trait1: I; trait2: C"
prop_better2$catg[prop_better2$catg == "10"] = "trait1: C; trait2: I"
prop_better2$catg[prop_better2$catg == "11"] = "trait1: C; trait2: C"
prop_better2$catg = factor(prop_better2$catg,
                           levels = c("trait1: I; trait2: I",
                                      "trait1: C; trait2: I",
                                      "trait1: I; trait2: C",
                                      "trait1: C; trait2: C"))
prop_better2$Models[prop_better2$Models == "PLMM"] = "PLMM better"
prop_better2$Models[prop_better2$Models == "LMM"] = "LMM better"

# data to plot power
type_errs2 = sim.FE2 %>% 
  group_by(trt, beta5, catg) %>%
  summarise(nitem = max(item),
            mid.est = median(estimate),
            ave.est = mean(estimate),
            ave.est.low = mean(estimate) - 1 * sd(estimate),
            ave.est.high = mean(estimate) + 1 * sd(estimate),
            reject.prop = sum(Pvalue <= 0.05)/nitem) %>%
  ungroup() 
type_errs2$catg[type_errs2$catg == "00"] = "trait1: I; trait2: I"
type_errs2$catg[type_errs2$catg == "01"] = "trait1: I; trait2: C"
type_errs2$catg[type_errs2$catg == "10"] = "trait1: C; trait2: I"
type_errs2$catg[type_errs2$catg == "11"] = "trait1: C; trait2: C"
type_errs2$catg = factor(type_errs2$catg,
                         levels = c("trait1: I; trait2: I",
                                    "trait1: C; trait2: I",
                                    "trait1: I; trait2: C",
                                    "trait1: C; trait2: C"))
type_errs2 = mutate(type_errs2, models = trt)
type_errs2$models[type_errs2$models == "no_phylo"] = "LMM"
type_errs2$models[type_errs2$models == "phylo"] = "PLMM"

# function to plot fig. 4 ----
est_se_plt = function(nsp = c(50, 20, 40, 60, 80), nitem = NULL, alpha = 1, 
                      trim = FALSE, trait2_C = TRUE, text1_x = 0.765, text2_x = 0.75,
                      text1_y = 0.07, text2_y = 0.05){
  nspp_est_se_0 = filter(sim.FE1, nspp == nsp)
  if(!is.null(nitem)){ # to plot less points
    idex = sample(1:max(nspp_est_se_0$item), size = nitem, replace = F)
    nspp_est_se = filter(nspp_est_se_0, item %in% idex)
  } 
  nspp_est_se$catg[nspp_est_se$catg == "00"] = "trait1: I; trait2: I"
  nspp_est_se$catg[nspp_est_se$catg == "01"] = "trait1: I; trait2: C"
  nspp_est_se$catg[nspp_est_se$catg == "10"] = "trait1: C; trait2: I"
  nspp_est_se$catg[nspp_est_se$catg == "11"] = "trait1: C; trait2: C"
  nspp_est_se$catg = factor(nspp_est_se$catg,
                            levels = c("trait1: I; trait2: I",
                                       "trait1: C; trait2: I",
                                       "trait1: I; trait2: C",
                                       "trait1: C; trait2: C"))
  nspp_est_se$trt[nspp_est_se$trt == "phylo"] = "PLMM"
  nspp_est_se$trt[nspp_est_se$trt == "no_phylo"] = "LMM"
  nspp_est_se = mutate(nspp_est_se, Significant = Pvalue < 0.05, Model = trt) 
  
  nspp_est_se_err = group_by(mutate(nspp_est_se_0, Model = trt), catg, Model) %>% 
    summarise(err = round(mean(Pvalue < 0.05), 2)) %>% 
    spread(key = Model, value = err) %>% 
    mutate(text1 = paste0("LMM = ", no_phylo), text2 = paste0("PLMM = ", phylo)) %>% 
    ungroup() %>% select(catg, text1, text2)
  
  nspp_est_se_err$catg[nspp_est_se_err$catg == "00"] = "trait1: I; trait2: I"
  nspp_est_se_err$catg[nspp_est_se_err$catg == "01"] = "trait1: I; trait2: C"
  nspp_est_se_err$catg[nspp_est_se_err$catg == "10"] = "trait1: C; trait2: I"
  nspp_est_se_err$catg[nspp_est_se_err$catg == "11"] = "trait1: C; trait2: C"
  nspp_est_se_err$catg = factor(nspp_est_se_err$catg,
                            levels = c("trait1: I; trait2: I",
                                       "trait1: C; trait2: I",
                                       "trait1: I; trait2: C",
                                       "trait1: C; trait2: C"))
  
  if(trim){
    if(trait2_C){
      nspp_est_se = filter(nspp_est_se, catg %in% c("trait1: I; trait2: C", 
                                                    "trait1: C; trait2: C")) %>% 
        mutate(catg = droplevels(catg))
      nspp_est_se_err = filter(nspp_est_se_err, catg %in% c("trait1: I; trait2: C", 
                                                            "trait1: C; trait2: C"))
    } else {
      nspp_est_se = filter(nspp_est_se, catg %in% c("trait1: I; trait2: I", 
                                                    "trait1: C; trait2: I")) %>% 
        mutate(catg = droplevels(catg))
      nspp_est_se_err = filter(nspp_est_se_err, catg %in% c("trait1: I; trait2: I", 
                                                            "trait1: C; trait2: I"))
    }
  }
  
  plt = ggplot(nspp_est_se, aes(x = abs(estimate), y = SE)) +
    geom_point(aes(color = Model), size = 1.5, shape = 1) +
    geom_abline(intercept = 0, slope = 1/1.96, linetype = 2) +
    geom_point(data = filter(nspp_est_se, Significant == TRUE), 
               aes(color = Model), alpha = alpha, size = 1.5, shape = 19) +
    facet_wrap(~catg, ncol = 2) +
    geom_text(data = nspp_est_se_err, aes(x = text1_x, y = text1_y, label = text1), 
              colour = "darkgoldenrod1", inherit.aes = FALSE, parse = FALSE) +
    geom_text(data = nspp_est_se_err, aes(x = text2_x, y = text2_y, label = text2), 
              colour = "blue", inherit.aes = FALSE, parse = FALSE) +
    labs(#title = paste0("nspp = ", nsp, "; true beta5 = 0"),
         x = expression(paste("Absolute value of estimates of ", beta[5], 
                              " (trait-environment interaction term)")),
         y = "Estimated standard error") + 
    # background_grid(major = "xy", minor = "none") +
    scale_color_manual(values = c("darkgoldenrod1", "blue")) +
    theme(legend.position = c(0.6, 0.9))
  print(plt)
}


## Fig 1: proportion of better estimations PLMM vs LMM ----
better_1_plot = filter(prop_better1, Models == "PLMM better") %>% 
  ggplot(aes(x = nspp, y = prop, color = catg, group = catg)) +
  geom_point(size = 2) + geom_line() + ylim(0, 1) +
  labs(y = "Relative frequency that PLMM is better", 
       x = "Number of species") +
  theme(legend.title = element_blank(), legend.position = c(0.85, 0.20)) +
  background_grid(major = "y", minor = "none")
  # 
  # ggplot(prop_better1, aes(x = nspp, y = prop)) +
  # geom_bar(stat = "identity", aes(fill = Models)) +
  # geom_hline(yintercept = 0.5, linetype = 2) +
  # facet_wrap(~catg, ncol = 1) +
  # scale_fill_manual(values = c("darkgoldenrod1", "blue")) +
  # labs(y = "Relative frequency", 
  #      x = "Number of species") +
  # theme(legend.title = element_blank(), legend.position = "top")
better_2_plot = filter(prop_better2, Models == "PLMM better") %>% 
  ggplot(aes(x = real.val, y = prop, color = catg, group = catg)) +
  geom_point(size = 2) + geom_line() + ylim(0, 1) +
  labs(x = expression(paste("True value of ", beta[5], 
                            " in Eq. 1 (trait-environment interaction term)")),
       y = "") +
  theme(legend.position = "none")+
  background_grid(major = "y", minor = "none")
  # 
  # ggplot(prop_better2, aes(x = real.val, y = prop)) +
  # geom_bar(stat = "identity", aes(fill = Models)) +
  # geom_hline(yintercept = 0.5, linetype = 2) +
  # facet_wrap(~catg, ncol = 1) +
  # scale_fill_manual(values = c("darkgoldenrod1", "blue")) +
  # labs(x = "True value of beta 5 in eq. 1",
  #      y = "") +
  # theme(legend.title = element_blank(), legend.position = "top")
plot_grid(better_1_plot, better_2_plot, labels = c("(A)", "(B)"))
if(!dir.exists("../figures")) dir.create("../figures")
ggsave(filename = "../figures/better_v2.pdf", width = 12, height = 6)

## Fig 2: biases of PLMM vs LMM ----
fig2_p = ggplot(mutate(type_errs, nspp = as.character(nspp)),
       aes(x = nspp, y = ave.est, group = models, color = models)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = ave.est.low, ymax = ave.est.high),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~catg) +
  theme(legend.position = "top", legend.title = element_blank()) +
  background_grid(major = "xy", minor = "none") +
  labs(x = "Number of species", 
       y = "Mean ± SD of the estimates") +
  scale_color_manual(values = c("darkgoldenrod1", "blue"))
dat = data.frame(x = rep(0.7, 4), y = rep(0.4, 4), catg = factor(unique(type_errs$catg),
                                                                 levels = c("trait1: I; trait2: I",
                                                                            "trait1: C; trait2: I",
                                                                            "trait1: I; trait2: C",
                                                                            "trait1: C; trait2: C")),
                 labs = c("(A)", "(C)", "(B)", "(D)"), models = "LMM")
fig2_p + geom_text(aes(x, y, label = labs, group = NULL), data = dat, color = "black", fontface = "bold")
ggsave(filename = "../figures/biases_1_2.pdf", width = 9, height = 8.1)

## Fig. 3: type I and power ----
type_1_plot = ggplot(type_errs,
       aes(x = nspp, y = reject.prop, color = models)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2) +
  facet_grid(catg~.) +
  theme(legend.position = c(0.8, 0.95),
        legend.title = element_blank(),
        strip.background = element_blank(),
        # panel.border = element_rect(color = "white"),
        strip.text.y = element_text(color = "white")) +
  background_grid(major = "xy", minor = "none") +
  labs(x = "Number of species",
       y = "Proportion of significant results") +
  scale_color_manual(values = c("darkgoldenrod1", "blue"))

type_2_plot = ggplot(type_errs2,
       aes(x = beta5, y = reject.prop, color = trt)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2) +
  facet_grid(catg ~ .) +
  theme(legend.position = "null") +
  labs(x = expression(paste("True value of ", beta[5], 
                            " in Eq. 1")),
       y = "") + 
  background_grid(major = "xy", minor = "none")+
  scale_color_manual(values = c("darkgoldenrod1", "blue"))

type_power_plot = plot_grid(type_1_plot, type_2_plot, labels = c("(A)", "(B)"))
type_power_plot
ggsave(filename = "../figures/type_power.pdf", width = 7, height = 10)

## Fig. 4 se vs abs(est) ----
est_se_plt(50, 300, alpha = 0.8, trim = T, text1_x = 0.763, text2_x = 0.75)
# this figure may be slighlty different from our paper as it randomly selected 300 points out of 1000.
ggsave(filename = "../figures/se_est_sp50_sample_300_2.pdf", height = 5, width = 10)

## Fig. A1: biases of PLMM vs LMM ----
ggplot(type_errs2,
       aes(x = models, y = ave.est, color = models)) +
  geom_hline(aes(yintercept = beta5), linetype = 2) +
  geom_pointrange(aes(ymin = ave.est.low, ymax = ave.est.high)) +
  facet_grid(beta5~catg, scales = "free") +
  labs(x = "Models",
       y = "Mean ± SD of the estimates") +
  theme(legend.position = "null") +
  background_grid(major = "xy", minor = "none") +
  scale_color_manual(values = c("darkgoldenrod1", "blue"))
ggsave(filename = "../figures/biases_2_2.pdf", width = 7, height = 8)

## Fig. A2 se vs abs(est) ----
est_se_plt(50, 300, alpha = 0.8, trim = T, trait2_C = F, text1_x = 0.386, text2_x = 0.38, text1_y = 0.06)
# this figure may be slighlty different from our paper as it randomly selected 300 points out of 1000.
ggsave(filename = "../figures/se_est_sp50_sample_300_appendix_2.pdf", height = 5, width = 10)
