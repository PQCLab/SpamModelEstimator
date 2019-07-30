% SPAM Model Estimation
[data, opt] = sme_init('example_ninf_mg.mat',struct());
opt.Fref = sme_loglik(data.spamtrue,data);
result = sme_solve(data, opt);