function [edges] = getNevTriggerEdges(nev, varargin)
% GETNEVTRIGGEREDGES Returns trigger UP and DOWN events as a Nx2
%   double matrix (EDGES).
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

options = struct(...
    'firstedgeup',true,...
    'window',[]);
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

events = double(nev.Data.SerialDigitalIO.TimeStamp);
amps = double(nev.Data.SerialDigitalIO.UnparsedData);

% clean up amplitude and remove any doubles
upthresh = mean(amps);
amps(amps > upthresh) = max(amps);
amps(amps <= upthresh) = min(amps);
deltaidx = find(abs(diff(amps))<1e-3);
if any(deltaidx)
    amps(deltaidx+1) = [];
    events(deltaidx+1) = [];
end

% restrict to a time window
if ~isempty(options.window)
    if ~isnumeric(options.window) || numel(options.window) ~= 2
        error('Sample limit must be a 2-element double.');
    end
    tmask = events >= options.window(1) & events <= options.window(2);
    events = events(tmask);
    amps = amps(tmask);
end

edgesup = events(amps > upthresh);
edgesdn = events(amps <= upthresh);

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
