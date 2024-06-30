twoq_ll_draws <- get_draws("two_q-selected",model_out_dir="~/Downloads/orig_s3_fits/",vars=c("log_lik"))

r_eff_twoq <- relative_eff(exp(twoq_ll_draws),cores = getOption("mc.cores", 1)) #calculate relative ESS, which allows for better estimates of the PSIS effective sample sizes and Monte Carlo error
two_q_loo <- loo(twoq_ll_draws, r_eff = r_eff_twoq)

twoqsb_ll_draws <- get_draws("two_q_sidebias-selected",model_out_dir="~/Downloads/orig_s3_fits/",vars=c("log_lik"))

r_eff_twoqsb <- relative_eff(exp(twoqsb_ll_draws),cores = getOption("mc.cores", 1)) #calculate relative ESS, which allows for better estimates of the PSIS effective sample sizes and Monte Carlo error
two_q_sidebias_loo <- loo(twoqsb_ll_draws, r_eff = r_eff_twoqsb)

loo_compare(two_q_sidebias_loo,two_q_loo)


oneq_ll_draws <- get_draws("one_q-selected",model_out_dir="~/Downloads/orig_s3_fits/",vars=c("log_lik"))

r_eff_oneq <- relative_eff(exp(oneq_ll_draws),cores = getOption("mc.cores", 1)) #calculate relative ESS, which allows for better estimates of the PSIS effective sample sizes and Monte Carlo error
one_q_loo <- loo(oneq_ll_draws, r_eff = r_eff_oneq)

loo_compare(two_q_sidebias_loo,two_q_loo,one_q_loo)



