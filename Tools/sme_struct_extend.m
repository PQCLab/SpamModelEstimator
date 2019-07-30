function c = sme_struct_extend(a, b)

c = a;
f = fieldnames(b);
for i = 1:length(f)
    c.(f{i}) = b.(f{i});
end

end

