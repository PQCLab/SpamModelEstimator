function M = sme_chi_povm(input, povm, schemes)

m = size(schemes,1);
M = cell(1,m);
for i = 1:m
    scheme = schemes(i,:);
    for j = 1:size(povm{scheme(2)})
        M{i}(:,:,j) = kron(conj(input{scheme(1)}), povm{scheme(2)}(:,:,j));
    end
end

end

