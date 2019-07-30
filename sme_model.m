function model = sme_model(init, readout, gates, prep, meas, G_EG)
%SME_MODEL Generates SPAM-model (input states and POVM-operatora for tomography)
%INPUT:
%   init -      initialization state (density matrix)
%   readout -   readout POVM-operators (3D array, each slice in 3-rd dimentions is POVM operator
%   gates -     gates list (cell array, gates{j}.G - evolution matrix for j-th gate)
%   prep -      preparation gate sequences
%   meas -      measurement gate sequances
%   G_EG -      empty gate evolution matrix (optional)
%OUTPUT:
%   model -     SPAM-model

if nargin < 6
    G_EG = 1;
end

d = size(init,1);
model.input = cell(1,length(prep));
for i = 1:length(prep)
    model.input{i} = reshape(seq2gate(prep{i},gates)*init(:),d,d);
end
readout = [reshape(readout(:,:,1),[],1), reshape(readout(:,:,2),[],1)]';
model.povm = cell(1,length(meas));
for i = 1:length(meas)
    povm = readout*seq2gate(meas{i},gates)*G_EG;
    model.povm{i} = cat(3, reshape(povm(1,:)',d,d), reshape(povm(2,:)',d,d));
end

end

function G = seq2gate(seq, gates)
    G = 1;
    for j = 1:length(seq)
        G = gates{seq(j)}.G * G;
    end
end
