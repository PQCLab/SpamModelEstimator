function G = sme_kraus2evol(E)

G = 0;
for k = 1:size(E,3)
    G = G + kron(conj(E(:,:,k)), E(:,:,k));
end

end

