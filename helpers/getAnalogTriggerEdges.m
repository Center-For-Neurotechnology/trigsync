function [edges] = getAnalogTriggerEdges(d, varargin)
% GETANALOGTRIGGEREDGES Returns trigger UP and DOWN events as a Nx2
%   double matrix (EDGES) within the K-element double vector D.
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

options = struct(...
    'maxvalpc',0.8,...
    'upthreshpc',0.5,...
    'firstedgeup',true);
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

% cap analog value
d = single(d);
maxval = options.maxvalpc*max(d);
upthresh = options.upthreshpc*maxval;

% clean up analog value
d(d > maxval) = maxval;
d(d < 0.5*maxval) = 0;

% capture the down edges
flipped = upthresh-d;
flipped(end) = 0;
[~,edgesup] = findpeaks(d,'minpeakheight',0.7*upthresh);
[~,edgesdn] = findpeaks(flipped,'minpeakheight',0.7*upthresh);

% edge filtering
if options.firstedgeup
    master.edges = edgesup; slave.edges = edgesdn;
else
    master.edges = edgesdn; slave.edges = edgesup;
end
master.num = numel(master.edges);
slave.num = numel(slave.edges);

% make sure the first edge is the master
if slave.edges(1) < master.edges(1)
    slave.edges(1) = []; slave.num = numel(slave.edges);
end

if slave.num > master.num || master.num-slave.num > 1
    error('Abnormal trigger events detected (are there enough triggers?).');
end

% ensure that master/slave edges come in pairs
if master.num > slave.num
    master.edges = master.edges(1:slave.num); master.num = numel(master.edges);
end
edges = [master.edges; slave.edges]';
