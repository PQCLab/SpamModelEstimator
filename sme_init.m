function [data, opt] = sme_init(fname, options)
%GST_INIT Initialize GST solver with specified options
%   fname -     data filenmae (must bet put in Data directory)
%   options -   options list

if ~exist('rt_chi_reconstruct', 'file')
    error('The current package depends on PQCLab/RootTomo. Add the path to this library to the searching path.');
end
if ~exist('rm_init', 'file')
    error('The current package depends on PQCLab/RandomMutations. Add the path to this library to the searching path.');
end
if ~exist('sme_struct_extend', 'file')
    addpath('Tools');
end

% === OPTIONS ===
if nargin < 2
    options = struct();
end
opt = sme_struct_extend(struct( ...
    'Nmax',          1e6, ...         % Maximum number of iteration
    'tolX',          1e-9, ...        % Tolerance for parameters
    'tolF',          1e-10, ...       % Tolerance for objective function
    'FminsTol',      1e-8, ...        % Tolerance for using of fminsearch
    'FminsFreq',     100, ...         % Number of steps between fnimsearch
    'Display',       3, ...           % Display mode (0 - silent, 1 - steps, 2 - iterations, 3 - all)
    'DispIterLines', 3, ...           % Number of lines in iterations
    'Fref',          'auto', ...      % Reference loglik
    'ForceRestart',  false, ...       % Restart calculation
    'irest',         true, ...        % Estimate initialization and readout
    'gest',          'all', ...       % Gates to estimate (array, all - estimate all tomo gates except empty)
    'rank',          1:4 ...          % Testing gates ranks (array)
), options);

opt.d = 2;
opt.d2 = 4;
opt.rank = sort(opt.rank);
opt.fdata = strcat('Data/', fname);
opt.fresult = strcat('Results/', fname);
if ~isfile(opt.fdata)
    error('Error. File %s not found.', opt.fdata);
end
load(opt.fdata, 'data');
if ~exist('Results', 'dir')
	mkdir('Results');
end

% === PREPARE RESULT ===
if (isfield(opt,'ForceRestart') && opt.ForceRestart) || ~exist(opt.fresult, 'file')
    result.status = 0;
    result.time = [];
    result.rank = opt.rank;
    result.init = data.init;
    result.readout = data.readout;
    [result.init_v, result.readout_v] = sme_irvec(result.init, result.readout);
    result.gates = data.gates;
    for it = 1:length(result.gates)
        result.gates{it}.G = sme_kraus2evol(result.gates{it}.U);
        result.gates{it}.chi = sme_kraus2chi(result.gates{it}.U);
    end
    result.prep = data.prep;
    result.meas = data.meas;
    result.model = sme_model(result.init, result.readout, result.gates, result.prep, result.meas);
    resultr = cell(1,length(opt.rank));
    save(opt.fresult, 'result', 'resultr');
else
    load(opt.fresult, 'result');
    if length(opt.rank) ~= length(result.rank) || norm(opt.rank - result.rank) > 0
        error('You have changed rank. Set option ForceRestart=true to restart calculations');
    end
end

% === PREPARE DATA ===
opt.tomoempty = 0;
opt.tomoest = [];
data.circuits = {};
data.schemes = [];
data.clicks = [];
data.nshots = [];
for it = 1:length(data.tomo)
    tomo = data.tomo{it};
    gate = data.gates{tomo.gate};
    if isfield(gate, 'empty') && gate.empty
        opt.tomoempty = it;
    elseif (ischar(opt.gest) && strcmp(opt.gest, 'all')) || ~isempty(find(opt.gest == tomo.gate, 1))
        opt.tomoest = [opt.tomoest, it];
    end
    tomo.nshots = sum(tomo.clicks, 2);
    ind = find(tomo.nshots == 0);
    if ~isempty(ind)
        warning('Empty results ignored: gate %d, schemes %d', tomo.gate, ind);
        tomo.schemes(ind,:) = [];
        tomo.clicks(ind,:) = [];
        tomo.nshots(ind) = [];
    end
    data.schemes = [data.schemes; tomo.schemes];
    data.clicks = [data.clicks; tomo.clicks];
    data.nshots = [data.nshots; tomo.nshots];
    
    schemes = mat2cell(tomo.schemes, ones(1,length(tomo.schemes)));
    circuits = cellfun(@(scheme) [data.prep{scheme(1)}, tomo.gate, data.prep{scheme(2)}], schemes, 'UniformOutput', false);
    data.circuits = [data.circuits; circuits];
    
    data.tomo{it} = tomo;
end
if opt.tomoempty == 0
    error('You must provide an empty gate results');
end
etgate = data.tomo{opt.tomoempty}.gate;
data.circuits = cellfun(@(circuit) circuit(circuit ~= etgate), data.circuits, 'UniformOutput', false);
opt.gateest = cellfun(@(tomo) tomo.gate, data.tomo(opt.tomoest));
opt.L = sum(data.nshots)*length(opt.tomoest);
if strcmp(opt.Fref, 'auto')
    opt.Fref = sme_loglik(result,data);
end

end

