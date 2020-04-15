function [x,y] = makeWaveformFromEdges(edges, varargin)
% MAKEWAVEFORMFROMEDGES
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

options = struct(...
    'presamples',10,...
    'firstedgeup',true,...
    'upval',1,...
    'fs',[]);
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

edges = edges';
edges = edges(:);
numpoints = 2*numel(edges);

x = nan(numpoints,1);
y = nan(numpoints,1);

nextval = options.firstedgeup;
for ii = 1:numpoints
    pidx = 1+(ii-1)*2;
    cidx = pidx+1; 
    
    if ii > 1, nextval = 1-lastval; end
    
    if pidx <= numpoints
        x(pidx) = edges(ii)-1; y(pidx) = 1-nextval;
    end
    if cidx <= numpoints
        x(cidx) = edges(ii); y(cidx) = nextval;
    end
    lastval = nextval;
end

% add some blank samples at the beginning
if ~isempty(options.presamples)
    x = [x(1)-options.presamples; x];
    y = [~y(1); y];
end

y = y*options.upval;
x = x-x(1);
% convert to time
if ~isempty(options.fs), x = x/options.fs; end
