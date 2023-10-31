tl_dat <- read.csv("/Users/dp/projects/ARL_bv/analysis_data/trial_level_data_all_subs_2023-04-20_19_56_01.csv")
tl_dat_aff <- tl_dat %>% filter(!is.na(valrat_z) & show_fres == 1)
summary(lm(valrat_z ~ out + block_val,tl_dat_aff))
cor(tl_dat_aff$valrat_z,tl_dat_aff$block_val)
cor(tl_dat_aff$valrat_z,tl_dat_aff$out)
tl_dat_aff_stay <- tl_dat_aff %>% filter(!is.na(stay))
cor.test(tl_dat_aff_stay$stay,tl_dat_aff_stay$block_val)
final_tl <- tl_dat_aff_stay %>% filter(id %in% c(79993,
                                                 82774,
                                                 87358,
                                                 80944,
                                                 87502,
                                                 89626,
                                                 88198))
cor.test(final_tl$stay,final_tl$block_val)
cor(final_tl$valrat_z,final_tl$block_val)

summary(glm(stay ~ out + valrat_z,tl_dat,family = "binomial"))

summary(glm(stay ~ out + valrat_z,tl_dat_aff,family = "binomial"))

hi_bv <- tl_dat %>% filter(id %in% c(89692,
                               87664,
                               87565,
                               88198,
                               79993))
hi_bv_fres <- hi_bv %>% filter(show_fres == 1)
summary(glm(stay ~ out + block_val,hi_bv_fres,family = "binomial"))
cor.test(hi_bv_fres$valrat_z,hi_bv_fres$block_val,na.rm=TRUE)
cor.test(hi_bv_fres$valrat_z,hi_bv_fres$out,na.rm=TRUE)


showfres <- tl_dat %>% filter(show_fres == 1)
summary(glm(stay ~ out + valrat_z,showfres,family = "binomial"))
nofres <- tl_dat %>% filter(show_fres == 0)
summary(glm(stay ~ out + valrat_z,nofres,family = "binomial"))
