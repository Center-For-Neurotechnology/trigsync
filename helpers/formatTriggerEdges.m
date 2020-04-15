function [edges] = formatTriggerEdges(edges, fs, trigset, varargin)
% FORMATTRIGGEREDGES Cleans up trigger edges to ensure consistency and
%	reliability (e.g. certain pulse edges appear first, depending on the
%	trigger type.
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

options = struct(...
    'autocorrect',true);
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

if strcmp(trigset,'presentation')
    % three pulse group (durations 0.1s, 0.2s, 0.3s), fixed amplitude
    %   * P1 and P2 separation encodes the hour
    %   * P2 and P3 separation encodes the minute
    %   * P3 to next group duration is fixed
    
    % arrange for the 100 ms pulse to be leading
    pulsems = round(1e3*diff(edges')/fs);
    idxlead = find(pulsems >= 98 & pulsems <= 102);
    
    % first, attempt to align the pulses by the first 100ms pulse
    % if there are no 100ms pulses, attempt to align by the 300ms pulse
    % TODO: handle this possibility!
    if isempty(idxlead)
        error('100 ms leading pulses not found for trigset %s.',trigset);
    end
    if idxlead(1) > 1
        edges = edges(idxlead(1):end,:);
        idxlead = idxlead - idxlead(1) + 1;
    end
    
    % when the time hits o'clock, the HH value will be zero (no pulse!)
    % so we need to fake this missing pulse for every leading pulse
    % this gets complicated if we are unlucky enough to go from
    % 23:59 to midnight, but hopefully this is sorted out
    idxmid = find(pulsems >= 198 & pulsems <= 202);
    nummidpulsesamples = round(200e-3*fs);
    numextrasamples = 4;

    if numel(idxlead) ~= numel(idxmid)
        nextslot = 1;
        numedges = size(edges,1);
        newedges = nan(3*numel(idxlead),2);
        for ii = 1:numedges
            newedges(nextslot,:) = edges(ii,:);
            nextslot = nextslot + 1;
            if ii >= numedges, break; end

            % add a second pulse straight after the leading pulse
            % if we don't get one where we expect to see one
            thisisleading = isAlmostEqual(diff(edges(ii,:))/fs,100e-3,'atol',1e-2);
            nextismiddle = isAlmostEqual(diff(edges(ii+1,:))/fs,200e-3,'atol',1e-2);
            if thisisleading && ~nextismiddle
                newedges(nextslot,1) = edges(ii,2)+numextrasamples;
                newedges(nextslot,2) = newedges(nextslot,1) + nummidpulsesamples;
                nextslot = nextslot + 1;
            end
        end
        edges = newedges;
    end

    % by this point, we expect there to always be three pulses in a
    % sequence: 100, 200, 300 ms; if not, correct using padding
    if options.autocorrect
        tempedges = edges';
        tempedges = tempedges(:);
        
        % ensure a multiple of 6 edges (3 pulses per group)
        if mod(numel(tempedges),6)
            tempedges = [tempedges; nan(6-mod(numel(tempedges),6),1)];
        end
        numgroups = numel(tempedges)/6;

        % go through each pulse group and correct errors
        for ii = 1:numgroups-1
            idx = (ii-1)*6 + (1:7);
            seps = diff(tempedges(idx))/fs;
            pwidths = round(1e3*seps([1 3 5])');

            % in this case add some nan padding to properly split the first
            % and second pulses, to push the 30s separation back into place
            if isAlmostEqual(pwidths,[300 300 100]) && ~isAlmostEqual(seps(end),29.7,'atol',1)
                tempedges = [tempedges(1:idx(1)-1); [nan nan]'; tempedges(idx(1):end)];
            end
        end
        
        edges = reshape(tempedges,2,[])';
    end
else
    error('Unrecognized trigset (%s).',trigset);
end
