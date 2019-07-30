% Test results
load('Data/example_ninf_mg.mat', 'data');
load('Results/example_ninf_mg.mat', 'result');
load('Results/example_ninf_mg_gst.mat', 'result_gst');

n = 10.^[1,2,3,4,5,6];
Nexp = 500;
n_asymp = 1e3;
U = [1 1; 1 -1]/sqrt(2);
gamma = 0.9;

E = zeros(2,2,4);
E(:,:,1) = sqrt(1-3*gamma/4)*U;
E(:,:,2) = sqrt(gamma/4)*[0,1;1,0];
E(:,:,3) = sqrt(gamma/4)*[0,-1j;1j,0];
E(:,:,4) = sqrt(gamma/4)*[1,0;0,-1];
chi_true = sme_kraus2chi(E); % true process chi-matrix

% True SPAM-model
model_true = sme_model(data.spamtrue.init,data.spamtrue.readout,data.spamtrue.gates,data.prep,data.meas);
P0_true = cat(3,model_true.input{:});
M0_true = model_true.povm;
POVM_true = rt_chi_protocol(P0_true, M0_true);

% Reconstructed SPAM-model (GST-library)
model_gst = sme_model(result_gst.init,result_gst.readout,result_gst.gates,data.prep,data.meas);
POVM_gst = rt_chi_protocol(cat(3,model_gst.input{:}), model_gst.povm);

% Reconstructed SPAM-model (SME-library)
POVM_sme = rt_chi_protocol(cat(3,result.model.input{:}), result.model.povm);

% Performing experiments
[F_true, F_gst, F_sme] = deal(zeros(length(n),Nexp));
for i = 1:length(n)
    for j = 1:Nexp
        fprintf('n = %d | i_exp = %d\n', n(i), j);
        clicks = rt_chi_simulate(chi_true, P0_true, M0_true, n(i));
        [F_true(i,j), F_gst(i,j), F_sme(i,j)] = reconstruct(clicks,POVM_true,POVM_gst,POVM_sme,n(i),'auto',chi_true);
    end
end

% Error bars
figure;
set(gca,'FontSize',14);
hold on; grid on;
plot_err(n, F_true, 'True model', '-*');
plot_err(n, F_gst, 'Reconstructed model (GST-lib)', '-o');
plot_err(n, F_sme, 'Reconstructed model (SME-lib)', '-s');
xlabel('$$n$$', 'Interpreter', 'latex');
ylabel('$$1-F$$', 'Interpreter', 'latex');
set(gca, 'XScale','log');
set(gca, 'YScale','log');
legend('show');

%% Theory and experiment
disp('=== Theory (true) vs experiment (sme) ===');
figure;
indn = 6;
h = histogram(1-F_true(indn,:));
hold on;
xmax = h.BinLimits(2);
dx = xmax / 100;
x = 0:dx:xmax;
d = rt_chi_theory(chi_true,POVM_true,n(indn));
P = chi2pdf_general(x,d);
N_dist = P * h.BinWidth * Nexp;
plot(x, N_dist, 'LineWidth', 2);
F_mean_theory = sum(d);
F_mean_experimaent = mean(1-F_true(indn,:));
table(F_mean_theory, F_mean_experimaent)

%%
disp('=== Asymptotic results ===');
clicks = rt_chi_simulate(chi_true, cat(3,model_true.input{:}), model_true.povm, n_asymp, true);
[F_true_asymp, F_gst_asymp, F_sme_asymp] = reconstruct(clicks,POVM_true,POVM_gst,POVM_sme,n_asymp,'auto',chi_true);
table(F_true_asymp, F_gst_asymp, F_sme_asymp)

% Helpers
function [F_true, F_gst, F_sme] = reconstruct(clicks,POVM_true,POVM_gst,POVM_sme,n,r,chi_true)
    chi_rec = rt_chi_reconstruct(clicks,POVM_true,n,r,false);
    F_true = rt_fidelity(chi_true, chi_rec);
    chi_rec = rt_chi_reconstruct(clicks,POVM_gst,n,r,false);
    F_gst = rt_fidelity(chi_true, chi_rec);
    chi_rec = rt_chi_reconstruct(clicks,POVM_sme,n,r,false);
    F_sme = rt_fidelity(chi_true, chi_rec);
end
function plot_err(n, F, name, style)
    dF = 1-F;
    dFm = mean(dF,2);
    dFq = quantile(dF, [0.25, 0.75], 2);
    errorbar(n, dFm, dFm-dFq(:,1), dFq(:,2)-dFm, style, 'LineWidth', 1.5, 'MarkerSize', 10, 'DisplayName', name);
end