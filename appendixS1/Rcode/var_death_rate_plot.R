list.files("rds_var_bd/")
test1 = readRDS("rds_var_bd/sp20site30_birth1death0.1_11.rds")
str(test1)
head(test1$est.fixed.eff)
max(test1$est.fixed.eff$item)

test1 = bind_rows(
readRDS("rds_var_bd/sp20site30_birth1death0.1_11.rds")$est.fixed.eff %>% mutate(nspp = 20, death_rate = 0.1) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp40site30_birth1death0.1_11.rds")$est.fixed.eff %>% mutate(nspp = 40, death_rate = 0.1) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp60site30_birth1death0.1_11.rds")$est.fixed.eff %>% mutate(nspp = 60, death_rate = 0.1) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp20site30_birth1death0.3_11.rds")$est.fixed.eff %>% mutate(nspp = 20, death_rate = 0.3) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp40site30_birth1death0.3_11.rds")$est.fixed.eff %>% mutate(nspp = 40, death_rate = 0.3) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp60site30_birth1death0.3_11.rds")$est.fixed.eff %>% mutate(nspp = 60, death_rate = 0.3) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp20site30_birth1death0.5_11.rds")$est.fixed.eff %>% mutate(nspp = 20, death_rate = 0.5) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp40site30_birth1death0.5_11.rds")$est.fixed.eff %>% mutate(nspp = 40, death_rate = 0.5) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp60site30_birth1death0.5_11.rds")$est.fixed.eff %>% mutate(nspp = 60, death_rate = 0.5) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp20site30_birth1death0.7_11.rds")$est.fixed.eff %>% mutate(nspp = 20, death_rate = 0.7) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp40site30_birth1death0.7_11.rds")$est.fixed.eff %>% mutate(nspp = 40, death_rate = 0.7) %>% filter(terms == "envi1:trait1"),
readRDS("rds_var_bd/sp60site30_birth1death0.7_11.rds")$est.fixed.eff %>% mutate(nspp = 60, death_rate = 0.7) %>% filter(terms == "envi1:trait1")) %>% 
  tbl_df()

group_by(test1, nspp, death_rate, trt) %>% summarise(n_sig = mean(Pvalue < 0.05)) %>% 
  ungroup() %>% 
  mutate(trt = revalue(trt, c("no_phylo" = "LMM", "phylo" = "PLMM")),
         death_rate = as.factor(death_rate)) %>% 
  ggplot(aes(x = nspp, y = n_sig, color = trt, shape = death_rate)) +
  geom_point() + geom_line() + geom_hline(yintercept = 0.05, linetype = 2) +
  labs(x = "Number of species", y = "Type I error", title = "Type I errors of LMM and PLMM with varing death rate, number of simulations = 200, nsite = 30, birth rate = 1.")
ggsave(filename = "var_death_rate_type_I.pdf", height = 8, width = 10)
