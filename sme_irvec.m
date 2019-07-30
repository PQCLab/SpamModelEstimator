function [init,readout] = sme_irvec(init,readout)
%SME_IRVEC Vectorizes (or devectorizes) init and readout

if size(init,2) == 1 % devectorize
    d = sqrt(length(init));
    init = reshape(init, d, d);
    readout = cat(3, reshape(readout(1,:)', d, d), reshape(readout(2,:)', d, d));
else % vectorize
    init = reshape(init, [], 1);
    readout = [reshape(readout(:,:,1),[],1), reshape(readout(:,:,2),[],1)]';
end

end

