function [edges] = getAnalogTriggerEdges(d, varargin)
% GETANALOGTRIGGEREDGES Returns trigger UP and DOWN events as a Nx2
%   double matrix (EDGES) within the K-element double vector D.
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

options = struct(...
    'maxvalpc',0.8,...
    'upthreshpc',0.5,...
    'upval',1,...
    'trigset',[],...
    'fs',[],...
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

% cap and clean analog value
d = single(d);
upthresh = options.maxvalpc*max(d);
d(d < upthresh) = 0;
d(d >= upthresh) = options.upval;

% capture the down edges
flipped = options.upval-d;
flipped(end) = 0;
[~,edgesup] = findpeaks(d,'minpeakheight',0.5*options.upval);
[~,edgesdn] = findpeaks(flipped,'minpeakheight',0.5*options.upval);

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

% sometimes analog triggers can be encoded as edge (i.e. digital) events
% this manifests as every pulse being one sample long
if numel(unique(diff(edges'))) == 1
    if strcmp(options.trigset,'presentation')
        if isempty(options.fs)
            error('Sample rate (fs) required to handle this trigset case.');
        end
        
        % in this case, the events are UP edges
        % find the 30s delay (will either have an index of 1, 2 or 3 and from there determine the sequence start
        dt = diff(master.edges)/options.fs;
        longidx = find(abs(dt-30) < 0.1,1);
        if longidx < 3
            master.edges = master.edges(longidx+1:end);
        end
        
        % make sure every rise has its fall
        seqoffsets = (1+mod(0:(numel(master.edges)-1),3))*100e-3;
        seqoffsets = round(seqoffsets*options.fs);
        edges = [master.edges; master.edges+seqoffsets]';
    else
        error('Digital-like trigger edges detected, but no trigset parameter supplied.');
    end
end
