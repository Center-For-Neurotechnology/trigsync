function result = isAlmostEqual(A, B, varargin)
% ISALMOSTEQUAL Evaluates whether the elements of A and B are numerically 
%   close to each other. Based on numpy's ISCLOSE function.
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

options = struct(...
    'atol',1e-8,...
    'rtol',1e-5);
paramNames = fieldnames(options);

nArgs = length(varargin);
if round(nArgs/2) ~= nArgs/2
    error('Name/value input argument pairs required.')
end

% {name; value} pairs
for pair = reshape(varargin,2,[])
    thisParam = lower(pair{1});
    if any(strcmp(thisParam,paramNames))
        options.(thisParam) = pair{2};
    else
        error('%s is not a recognized parameter name.',thisParam)
    end
end

% -------------------------------------------------------------------------

result = false;

% check dimensions
if ndims(A) ~= ndims(B), return; end
for ii = 1:ndims(A)
    if size(A,ii) ~= size(B,ii), return; end
end

% now check the values
A = A(:);
B = B(:);
for ii = 1:numel(A)
    if abs(A(ii)-B(ii)) >= (options.atol + options.rtol*abs(B(ii))), return; end
end

result = true;
