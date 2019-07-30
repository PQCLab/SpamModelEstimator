function result = sme_solve(data, opt)
%SME_SOLVE Finding SPAM model
%   data -  Experimental data
%   opt -   Optimization options
%TODO conjugence error
load(opt.fresult, 'result', 'resultr');

rm_init;
rmparams.display = 0;
rmparams.p_max = 0; 
rmparams.n_stall = 20;

if opt.Display > 0
    h = sme_print('EGT SPAM-errors: initialization\n');
end
if result.status < 1
    tic;
    if opt.Display > 0
        h = sme_print('EGT SPAM-errors: calculation\n', h);
    end
    tomo = data.tomo{opt.tomoempty};
    M = sme_chi_povm(result.model.input, result.model.povm, tomo.schemes);
    chi = rt_chi_reconstruct(tomo.clicks', M, tomo.nshots', 'auto', false)*opt.d;
    result.model_egt = sme_model(result.init, result.readout, result.gates, result.prep, result.meas, sme_chi2evol(chi));
    result.status = 1;
    result.time = [result.time, toc];
    resultr(:) = {result};
    save(opt.fresult, 'result', 'resultr');
    if opt.Display > 0
        h = sme_print('EGT SPAM-errors: DONE\n', h);
    end
elseif opt.Display > 0
    h = sme_print('EGT SPAM-errors: LOADED\n', h);
end

% SME
options = optimset('Display','none','MaxFunEvals',3e4,'MaxIter',2e4,'TolX',1e-10,'TolFun',1e-10); %fminsearch
e0 = zeros(opt.d2,opt.d2,length(opt.tomoest)+1);
e0(1,1,end) = max(0.00001, result.init_v(4));
e0(2,1,end) = max(0.00001, result.readout_v(2,1));
e0(3,1,end) = max(0.00001, result.readout_v(1,4));
opt.circind = cell(1,length(opt.tomoest));
opt.circcont = cell(1,length(opt.tomoest));
for ie = 1:length(opt.tomoest)
    opt.circind{ie} = cellfun(@(circuit) find(circuit == data.tomo{opt.tomoest(ie)}.gate), data.circuits, 'UniformOutput', false);
    opt.circcont{ie} = ~cellfun('isempty', opt.circind{ie});
end
opt.circcont = cell2mat(opt.circcont);

for ir = 1:length(opt.rank)
    r = opt.rank(ir);
    if opt.Display > 0
        fprintf('=== Processing Rank-%d ===\n', r);
    end
    
    if opt.Display > 0
        h = sme_print('Zero approximation by EGT: initialization\n');
    end
    spam = resultr{ir};
    if spam.status < 2
        tic;
        spam.gateest = opt.gateest;
        spam.rest = r;
        if opt.Display > 0
            h = sme_print('Zero approximation by EGT: calculation\n', h);
        end
        for ie = 1:length(opt.tomoest)
            tomo = data.tomo{opt.tomoest(ie)};
            M = sme_chi_povm(result.model_egt.input, result.model_egt.povm, tomo.schemes);
            chi = rt_chi_reconstruct(tomo.clicks', M, tomo.nshots', r, false)*opt.d;
            spam.gates{tomo.gate}.chi = chi;
            spam.gates{tomo.gate}.G = sme_chi2evol(chi);
        end
        spam.status = 2;
        spam.time = [spam.time, toc];
        resultr{ir} = spam;
        save(opt.fresult, 'result', 'resultr');
        if opt.Display > 0
            h = sme_print('Zero approximation by EGT: DONE\n', h);
        end
    elseif opt.Display > 0
        h = sme_print('Zero approximation by EGT: LOADED\n', h);
    end
    
    if opt.Display > 0
        h = sme_print('Maximum likelihood: initialization\n');
    end
    if spam.status < 3
        tic;
        e = e0;
        for ie = 1:length(opt.tomoest)
            [U,D] = svd(spam.gates{data.tomo{opt.tomoest(ie)}.gate}.chi);
            d = diag(D); d((r+1):end) = 0; d = d / sum(d) * opt.d;
            e(:,:,ie) = U*diag(sqrt(d));
        end
        if opt.Display > 0
            h = sme_print('Maximum likelihood: calculation\n', h);
            h.lines = opt.DispIterLines;
        end
        x0 = e; x1 = e; z1 = e; t0 = 0; t1 = 1; alp = 1/opt.L; F_x1 = 0; spamX = spam;
        [opt.ind_re, opt.ind_im] = parind(opt,r);
        for i = 1:opt.Nmax
            if mod(i,opt.FminsFreq) == 0 && abs(dF) > opt.FminsTol
                if opt.Display > 1
                    h = sme_print(sprintf('%d | %d | Fminsearch...\n', r, i), h);
                end
                x = fminsearch(@(x) fminsearchFun(x, opt, spamX, data, r), e2x(x1,opt), options);
                x1 = x2e(x,opt);
                if ~opt.irest
                    x1(:,:,end) = e(:,:,end);
                end
                continue;
            end
            y1 = x1 + t0/t1*(z1-x1) + (t0-1)/t1*(x1-x0);
            spamY = adjust_spam(spamX,y1,data.tomo,opt);
            z2 = proj(y1 + alp * grad(y1, spamY, data, opt, r), opt, r);
            v2 = proj(x1 + alp * grad(x1, spamX, data, opt, r), opt, r);
            spamZ = adjust_spam(spamY,z2,data.tomo,opt);
            spamV = adjust_spam(spamX,v2,data.tomo,opt);
            F_z2 = sme_loglik(spamZ, data);
            F_v2 = sme_loglik(spamV, data);
            if F_z2 >= F_v2
                x2 = z2;
                F_x2 = F_z2;
                spamX = spamZ;
            else
                x2 = v2;
                F_x2 = F_v2;
                spamX = spamV;
            end
            x0 = x1; x1 = x2; z1 = z2;
            t0 = t1; t1 = (sqrt(4*t0^2+1)+1)/2;
            dx = norm(x1(:)-x0(:));
            dF = F_x2 - F_x1;
            F_x1 = F_x2;
            if opt.Display > 1
                h = sme_print(sprintf('%d | %d | dx = %.4e | dF = %.4e | dFt = %.4e |\n', r, i, dx, dF, opt.Fref-F_x1), h);
            end
            if (dx < opt.tolX && abs(dF) < opt.tolF)
                break;
            end
        end
        spam = spamX;
        [spam.init, spam.readout] = sme_irvec(spam.init_v, spam.readout_v);
        spam.time = [spam.time, toc];
        spam.status = 3;
        resultr{ir} = spam;
        save(opt.fresult, 'result', 'resultr');
        if opt.Display > 0
            h.lines = 1;
            h = sme_print('Maximum likelihood: DONE\n', h);
        end
    elseif opt.Display > 0
        h = sme_print('Maximum likelihood: LOADED\n', h);
    end
    
    if opt.Display > 0
        h = sme_print('Gauge optimization: initialization\n');
    end
    if spam.status < 4
        tic;
        if r > 1
            if opt.Display > 0
                h = sme_print('Gauge optimization: STEP 1/3 (this may take a while)\n', h);
            end
            rmparams.n_pop = 30;
            rmparams.n_des = 40;
            rmparams.eps = 1e-7;
            scale = 0.01;
            x = [0;0;0;1;0;0;0;0]; a0 = spam.init_v(4); a = spam.readout_v(2,1); b = spam.readout_v(1,4);
            for i = 1:10
                [f,x] = OptimizeRM(@(x) optimizermFun(x, spam.gates, spam.gateest, a0, a, b), 8, x, scale);
                if f < 100
                    break;
                end
                scale = scale / 10;
            end
            if opt.Display > 0
                h = sme_print('Gauge optimization: STEP 2/3 (this may take a while)\n', h);
            end
            rmparams.eps = 1e-8;
            [~,x] = OptimizeRM(@(x) optimizermFun(x, spam.gates, spam.gateest, a0, a, b), 8, x, scale/100);
            if opt.Display > 0
                h = sme_print('Gauge optimization: STEP 3/3 (this may take a while)\n', h);
            end
            rmparams.n_pop = 40; 
            rmparams.n_des = 60;
            rmparams.eps = 1e-9; 
            [~,x] = OptimizeRM(@(x) optimizermFun(x, spam.gates, spam.gateest, a0, a, b), 8, x, scale/100000);
            spam.Gauge = x2gauge(x,a0,a,b);
            spam.gates = gauge_gates(spam.gates, spam.gateest, spam.Gauge);
            spam.init_v = spam.Gauge * spam.init_v;
            spam.readout_v = spam.readout_v / spam.Gauge;
            [spam.init, spam.readout] = sme_irvec(spam.init_v, spam.readout_v);
            if opt.Display > 0
                h = sme_print('Gauge optimization: DONE\n', h);
            end
        elseif opt.Display > 0
        	h = sme_print('Gauge optimization: no need for rank-1\n', h);
        end
        spam.time = [spam.time, toc];
        spam.status = 4;
        resultr{ir} = spam;
        save(opt.fresult, 'result', 'resultr');
    elseif opt.Display > 0
        h = sme_print('Gauge optimization: LOADED\n', h);
    end
    
    [spam.adeq.pval, spam.adeq.chi2, spam.adeq.df] = sme_adeq(spam,data,opt,r);
    if opt.Display > 0
        fprintf('Stats: time = %d sec, p-value = %.4f, chi2 = %.4e, df = %d\n', round(sum(spam.time)), spam.adeq.pval, spam.adeq.chi2, spam.adeq.df);
    end
    
    spam.model = sme_model(spam.init, spam.readout, spam.gates, spam.prep, spam.meas);
    resultr{ir} = spam;
    result = best_result(resultr);
    save(opt.fresult, 'result', 'resultr');
end
end

function result = best_result(resultr)
    pvals = zeros(1,length(resultr));
    for i = 1:length(resultr)
        if isfield(resultr{i},'adeq')
            pvals(i) = resultr{i}.adeq.pval;
        end
    end
    indbest = find(pvals == max(pvals), 1);
    result = resultr{indbest};
    result.status = 5;
end
% === For proximal descend ===
function spam = adjust_spam(spam,e,tomo,opt)
    for ie = 1:length(opt.tomoest)
        chi = e(:,:,ie)*e(:,:,ie)';
        chi = chi / trace(chi) * opt.d;
        spam.gates{tomo{opt.tomoest(ie)}.gate}.chi = chi;
        spam.gates{tomo{opt.tomoest(ie)}.gate}.G = sme_chi2evol(chi);
    end
    if opt.irest
        spam.init_v = [1-e(1,1,end); 0; 0; e(1,1,end)];
        spam.readout_v = [1-e(2,1,end), 0, 0, e(3,1,end); e(2,1,end), 0, 0, 1-e(3,1,end)];
    end
end
function gr = grad(e,spam,data,opt,r)
    gr = zeros(size(e));
    gr_init = [-1; 0; 0; 1];
    gr_readout0 = [-1 0 0 0; 1 0 0 0];
    gr_readout1 = [0 0 0 1; 0 0 0 -1];
    gates = [spam.gates, struct('G', eye(4))];
    indlast = length(gates);
    for i = 1:length(data.circuits)
        Gc = sme_circuit2gate(data.circuits{i}, gates);
        Rm = Gc * spam.init_v;
        p = real(spam.readout_v * Rm)';
        coeff = data.clicks(i,:)./p - data.nshots(i);
        
        for ie = find(opt.circcont(i,:))
            for p = 1:opt.d2
                for q = 1:r
                    chi = zeros(opt.d2);
                    chi(:,p) = 2*e(:,q,ie);
                    gates{indlast}.G = sme_chi2evol(chi);
                    G = 0;
                    for ieg = opt.circind{ie}{i}
                        dcirc = data.circuits{i};
                        dcirc(ieg) = indlast;
                        G = G + sme_circuit2gate(dcirc, gates);
                    end
                    gr(p,q,ie) = gr(p,q,ie) + coeff*(spam.readout_v * G * spam.init_v);
                end
            end
        end
        
        if opt.irest
            gr(1,1,end) = gr(1,1,end) + coeff*real(spam.readout_v*Gc*gr_init);
            gr(2,1,end) = gr(2,1,end) + coeff*real(gr_readout0*Rm);
            gr(3,1,end) = gr(3,1,end) + coeff*real(gr_readout1*Rm);
        end
    end
end
function e = proj(e,opt,r)
    for ie = 1:length(opt.tomoest)
        ec = e(:,:,ie);
        Ue = zeros(opt.d2, opt.d);
        for i = 1:r
            Ue((1:opt.d)+(i-1)*opt.d, 1:opt.d) = reshape(ec(:,i), opt.d, opt.d);
        end
        [U,~,V] = svd(Ue,'econ');
        Ue = U*V';
        for i = 1:r
            E = Ue(((i-1)*opt.d+1):(i*opt.d),1:opt.d);
            ec(:,i) = reshape(E, opt.d2, 1);
        end
        e(:,:,ie) = ec;
    end
    
    if opt.irest
        e(1,1,end) = max(0.00001, min(0.4, real(e(1,1,end))));
        e(2,1,end) = max(0.00001, min(0.4, real(e(2,1,end))));
        e(3,1,end) = max(0.00001, min(0.4, real(e(3,1,end))));
    end
end
% === For fminsearch ===
function [ind_re, ind_im] = parind(opt, r)
    synd_e = zeros(opt.d2, opt.d2, length(opt.tomoest)+1);
    synd_e(:,1:r,1:(end-1)) = 1;
    ind_im = nonzeros(reshape(reshape((1:numel(synd_e))', size(synd_e)) .* synd_e, [], 1));
    if opt.irest
        synd_e(1:3,1,end) = 1;
    end
    ind_re = nonzeros(reshape(reshape((1:numel(synd_e))', size(synd_e)) .* synd_e, [], 1));
end
function x = e2x(e,opt)
    x = [real(e(opt.ind_re)); imag(e(opt.ind_im))];
end
function e = x2e(x,opt)
    e = zeros(opt.d2,opt.d2,length(opt.tomoest)+1);
    e(opt.ind_re) = e(opt.ind_re) + x(1:length(opt.ind_re));
    e(opt.ind_im) = e(opt.ind_im) + 1j*x((length(opt.ind_re)+1):end);
end
function F = fminsearchFun(x, opt, spam, data, r)
    e = x2e(x,opt);
    e = proj(e,opt,r);
    spam = adjust_spam(spam,e,data.tomo,opt);
    F = -sme_loglik(spam, data) + opt.Fref;
end
% === For gauge optimization ===
function Gauge = x2gauge(x,a0,a,b)
    w = x(1:2);
    p2max = min([1+(1-a0)/a0*w(1),  1-(1-b)/a*w(1), b/(1-a)*(1-w(1))]);
    w(2) = w(2) - (w(2)-p2max)*double(w(2) > p2max);
    A = zeros(2,4);
    A(1,1:3) = reshape(x(3:5),1,3) + 1j*reshape(x(6:8),1,3);
    A(2,1:3) = conj(A(1,[1,3,2]));
    A(:,end) = -(1-a0)/a0*A(:,1);
    Gauge = vertcat([1-w(1), 0, 0, w(2)], A, [w(1), 0, 0, 1-w(2)]);
end
function gates = gauge_gates(gates, gnum, Gauge)
    for ig = gnum
        gates{ig}.G = Gauge*gates{ig}.G/Gauge;
        gates{ig}.chi = sme_evol2chi(gates{ig}.G);
    end
end
function F = optimizermFun(x,gates,gnum,a0,a,b)
    gates = gauge_gates(gates, gnum, x2gauge(x,a0,a,b));
    F = 0; dn = 0;
    for ig = gnum
        d = eig(gates{ig}.chi);
        dn = dn + sum(d(d<0));
        F = F + norm(sme_kraus2evol(gates{ig}.U) - gates{ig}.G, 'fro');
    end
    F = F + 100*double(dn < 0) + 100*abs(dn) + (exp(100*abs(dn))-1);
end