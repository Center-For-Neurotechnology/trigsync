function [s] = computeTriggerSeparations(edges, fs, trigset, varargin)
% COMPUTETRIGGERSEPARATIONS(DATA, TTYPE) (wip)
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

options = struct(...
    'limit',[],...
    'infernan',true,...
    'roundunit',5e-2);
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

% transpose makes some operations easier
edges = edges'/fs;
numpulses = size(edges,2);
numedges = numel(edges);

% trigger type
if strcmp(trigset,'presentation')
    % three pulse group (durations 0.1s, 0.2s, 0.3s), fixed amplitude
    %   * P1 and P2 separation encodes the hour
    %   * P2 and P3 separation encodes the minute
    %   * P3 to next group duration is fixed
    numedgespergroup = 7;
    if numpulses < 4
        warning('Not enough pulses to compute %s separations.',trigset);
    end

    %numgroups = floor((numedges-1)/numsepspergroup);
    numgroups = floor(numedges/numedgespergroup);
else
    error('Unrecognized trigger type (%s).',trigset);
end

numsepspergroup = floor(numedgespergroup/2);
s = nan(numgroups,numsepspergroup);
for gg = 1:numgroups
    gidx = (gg-1)*(numedgespergroup-1) + (2:numedgespergroup);
    theseedges = diff(reshape(edges(gidx),2,[]));
    s(gg,:) = theseedges;
end

% round to requested units
if ~isempty(options.roundunit)
    s = round(s/options.roundunit)*options.roundunit;
end

nanidx = find(isnan(s(:)));
if options.infernan && ~isempty(nanidx)
    fprintf('\tInferring %d missing separations.\n',numel(nanidx));
    if strcmp(trigset,'presentation')
        % presentation triggers: seconds are constant, hours/min vary and
        % are inferred based on neighboring pulses
        for ii = 1:numel(nanidx)
            nidx = nanidx(ii);
            % avoid first/last groups
            if ~mod(nidx,numgroups) || ismember(nidx,[1 numel(s)])
                continue;
            end
            
            if nidx > 2*numgroups
                % seconds
                newval = 29.7;
            elseif nidx > numgroups
                % minutes (shouldn't happen!)
                isoclocked = diff(s([nidx-1,nidx]-numgroups)) < 0;
                if isoclocked
                    % TODO: check this!
                    newval = 1.0;
                else
                    % this is trickier; for now just take neighbor average
                    newval = mean(s([nidx-1,nidx+1]));
                end
            else
                % hours
                isoclocked = diff(s([nidx-1,nidx]+numgroups)) < 0;
                if isoclocked
                    % TODO: check this!
                    newval = mod(s(nidx-1)+0.05,0.6);
                else
                    newval = s(nidx-1);
                end
            end
            s(nidx) = newval;
        end
    end
end

% return no more than OPTIONS.LIMIT separations
if ~isempty(options.limit)
    limit = min([options.seplimit,numpulses]);
    s(limit+1:end,:) = [];
end
