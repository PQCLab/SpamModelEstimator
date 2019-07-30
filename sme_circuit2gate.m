function G = sme_circuit2gate(circuit, gates)

G = 1;
for ig = circuit
    G = gates{ig}.G * G;
end

end

