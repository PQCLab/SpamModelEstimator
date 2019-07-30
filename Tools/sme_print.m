function h = sme_print(str,h)

init = false;
if nargin < 2
    lines = 1;
    init = true;
end

if nargin == 2 && ~isstruct(h)
    lines = h;
    init = true;
end

if init
    h = struct('str', cell(1), 'lines', lines, 'nbytes', 0);
else
    fprintf(repmat('\b',1,h.nbytes));
end

h.str = [h.str, {str}];
if length(h.str) > h.lines
    h.str = h.str((end-h.lines+1):end);
end
h.nbytes = fprintf(strjoin(h.str,''));

end

