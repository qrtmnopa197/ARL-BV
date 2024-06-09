sign_flip <- function(vals){
  if(vals[1] < 0){
    vals[2] <- vals[2]*-1
  }
  return(vals[2])
}
pseudo_correct <- function(mod,model_out_dir,mu,sd){
  arr <- get_draws(mod,model_out_dir,c(mu,sd))
  return(apply(arr,c(1,2),sign_flip))
}


scd_nuis_fr_sens_mu <- pseudo_correct("temp_breakdown",model_out_dir,"nuis_fr_sens_mu","scd_nuis_fr_sens_mu")

quantile(scd_nuis_fr_sens_mu,c(.05,.5,.95))
mean(scd_nuis_fr_sens_mu < 0)


scd_ma_fr_sens_mu <- pseudo_correct("temp_breakdown",model_out_dir,"ma_fr_sens_mu","scd_ma_fr_sens_mu")

quantile(scd_ma_fr_sens_mu,c(.05,.5,.95))
mean(scd_ma_fr_sens_mu < 0)

scd_resid_fr_sens_mu <- pseudo_correct("temp_breakdown",model_out_dir,"resid_fr_sens_mu","scd_resid_fr_sens_mu")

quantile(scd_resid_fr_sens_mu,c(.05,.5,.95))
mean(scd_resid_fr_sens_mu < 0)

mean(scd_ma_fr_sens_mu > scd_rew_fr_sens_mu)
mean(scd_rew_fr_sens_mu > scd_resid_fr_sens_mu)


tb_scds <- get_draws("temp_breakdown",model_out_dir=model_out_dir,vars = c("scd_aff_fr_sens_mu","scd_rew_fr_sens_mu","scd_csens_mu"))
quantile(tb_scds[,,"scd_aff_fr_sens_mu"],c(.05,.5,.95))
quantile(tb_scds[,,"scd_rew_fr_sens_mu"],c(.05,.5,.95))
quantile(tb_scds[,,"scd_csens_mu"],c(.05,.5,.95))

spef_scds <- get_draws("temp_breakdown",model_out_dir=model_out_dir,vars = c("scd_ma_fr_sens_mu","scd_resid_fr_sens_mu","scd_nuis_fr_sens_mu"))

mean(tb_scds[,,"scd_aff_fr_sens_mu"] - tb_scds[,,"scd_rew_fr_sens_mu"] > 0)*100
mean(spef_scds[,,"scd_ma_fr_sens_mu"] - tb_scds[,,"scd_rew_fr_sens_mu"] > 0)*100

mean(tb_scds[,,"scd_rew_fr_sens_mu"] > 0)
mean(spef_scds[,,"scd_ma_fr_sens_mu"] > 0)
mean(spef_scds[,,"scd_resid_fr_sens_mu"] > 0)

quantile(spef_scds[,,"scd_ma_fr_sens_mu"],c(.05,.5,.95))
quantile(spef_scds[,,"scd_resid_fr_sens_mu"],c(.05,.5,.95))
quantile(spef_scds[,,"scd_nuis_fr_sens_mu"],c(.05,.5,.95))


two_q_scds <- get_draws("two_q",model_out_dir,vars=c("scd_rew_fr_sens_mu","scd_aff_fr_sens_mu"))
quantile(two_q_scds[,,"scd_rew_fr_sens_mu"],c(.05,.5,.95))
mean(two_q_scds[,,"scd_rew_fr_sens_mu"] > 0)
quantile(two_q_scds[,,"scd_aff_fr_sens_mu"],c(.05,.5,.95))
mean(two_q_scds[,,"scd_aff_fr_sens_mu"] > 0)
mean(two_q_scds[,,"scd_rew_fr_sens_mu"] > two_q_scds[,,"scd_aff_fr_sens_mu"])

scd_ma <- get_draws("breakdown",model_out_dir,vars=c("scd_ma_fr_sens_mu"))
quantile(scd_ma,c(0.05,.5,.95))

scd_rew <- get_draws("breakdown",model_out_dir,vars=c("scd_rew_fr_sens_mu"))
mean(scd_ma > scd_rew)
quantile(scd_rew,c(0.05,.5,.95))

scd_resid <- get_draws("breakdown",model_out_dir,vars=c("scd_resid_fr_sens_mu"))
quantile(scd_resid,c(0.05,.5,.95))
