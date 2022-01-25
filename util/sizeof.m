function nBytes = sizeof(precision)
%%% SIZEOF  return the number of bytes of a builtin data type.
error(nargchk(1, 1, nargin, 'struct'));
try
    d = regexp(str, '\d*','Match');
    nBits = str2num(d{1});
    nBytes = ceil(nBits / 8);
catch
    error('Unsupported class for finding size');
end
