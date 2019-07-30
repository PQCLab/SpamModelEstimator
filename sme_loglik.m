function logL = sme_loglik(spam,data)

if ~isfield(spam, 'init_v')
    [spam.init_v, spam.readout_v] = sme_irvec(spam.init, spam.readout);
end

logL = 0;
for i = 1:length(data.circuits)
    p = real(spam.readout_v * sme_circuit2gate(data.circuits{i}, spam.gates) * spam.init_v)';
    logp = log(p);
    logp(isnan(logp) | isinf(logp)) = 0;
    logL = logL + sum(data.clicks(i,:).*logp - p*data.nshots(i));
end

end
