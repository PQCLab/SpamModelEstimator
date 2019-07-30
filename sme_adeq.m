function [pval, chi2, df] = sme_adeq(spam,data,opt,r)

nu = ((2*opt.d2-r)*r - opt.d2)*length(opt.tomoest) + 3*double(opt.irest);
df = length(data.circuits) - nu;

if size(spam.init,2) > 1
    spam.init = spam.init(:);
    spam.readout = [reshape(spam.readout(:,:,1),[],1), reshape(spam.readout(:,:,2),[],1)]';
end

nE = zeros(size(data.clicks));
for i = 1:length(data.circuits)
    p = real(spam.readout * sme_circuit2gate(data.circuits{i}, spam.gates) * spam.init)';
    nE(i,:) = p*data.nshots(i);
end

nO = data.clicks;
chi2 = sum((nE-nO).^2./nE, 'all');
pval = vpa(gammainc(chi2/2,df/2,'upper'), 1000);

end

