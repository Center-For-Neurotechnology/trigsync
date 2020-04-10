function [t0] = locateCodedTrigger(targetpath, originpath, trigset, varargin)
% LOCATECODEDTRIGGER returns the sample of the first trigger in a TARGET
%   ephys file, and the corresponding sample of that first trigger within
%   an ORIGIN ephys file. These triggers are coded, meaning that their
%   value is related to their location in the recording.
%
%   TARGETPATH and ORIGINPATH can point to NSx, NEV, or EDF data, though
%   functionality has not yet been developed for some type combinations.
%   Please note that this is an alpha release! Feedback/testing is welcome.
%
%   For non-NEV data, a type of binary search is performed to locate TARGET
%   triggers within the ORIGIN structure. The ORIGIN search window
%   continually gets halved and shifted with each iteration.
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

options = struct(...
    'targetoffsettime',[],...
    'targettrigchan',[],...
    'origintrigchan',[],...
    'origintargetratio',2,...
    'nevpresamples',1e3,...
    'maxamppc',0.8,...
    'lower',[],...
    'upper',[],...
    'maxattempts',16,...
    'validate',false,...
    'searchlength',120);
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

% override default trigger channels if user-specified
targ.trigchan = options.targettrigchan;
orig.trigchan = options.origintrigchan;

% allow for direct NEV/NSx structure input
targ.nsx = []; targ.nev = []; targ.name = 'target'; targ.path = targetpath;
orig.nsx = []; orig.nev = []; orig.name = 'origin'; orig.path = originpath;

for trig = [targ,orig]
    if isstruct(trig.path) && isfield(trig.path,'MetaTags')
        ext = trig.path.MetaTags.FileExt;
        tpath = fullfile(trig.path.MetaTags.FilePath,...
            sprintf('%s%s',trig.path.MetaTags.Filename,ext));
        if strcmp(trig.name,'target')
            if isNsxPath(ext), targ.nsx = targetpath; else, targ.nev = targetpath; end
            targetpath = tpath;
            targ.path = tpath;
        else
            if isNsxPath(ext), orig.nsx = originpath; else, orig.nev = originpath; end
            originpath = tpath;
            orig.path = tpath;
        end
    end
end

% start sample for target trigger window
targ.t0 = 1;
targ.fs = getFileSampleRate(targetpath);
if ~isempty(options.targetoffsettime)
    targ.t0 = 1+floor(options.targetoffsettime*targ.fs);
end

% prepare target trigger window
targ.numsamples = round(options.searchlength*targ.fs);
targ.window = [targ.t0, targ.t0+targ.numsamples-1];

% get target trigger edges
if isNsxPath(targetpath)
    if isempty(targ.trigchan)
        targ.trigchan = getNsxTriggerChannel(targetpath,trigset);
    end
    if isempty(targ.nsx)
        targ.nsx = openNSx(targetpath,'channels',targ.trigchan,'duration',...
            targ.window(1):targ.window(2));
    end
    targ.data = targ.nsx.Data;
    targ.trigedges = getAnalogTriggerEdges(targ.data);
    targ.type = 'NSx';
elseif isNevPath(targetpath)
    if isempty(targ.nev)
        targ.nev = openNEV(targetpath,'nosave','nomat');
    end
    targ.trigedges = getNevTriggerEdges(targ.nev,'window',targ.window);
    targ.type = 'NEV';
elseif isEdfPath(targetpath)
    warning('This has not been tested!');
    targ.hdr = ft_read_header(targetpath);
    if isempty(targ.trigchan)
        targ.trigchan = getEdfTriggerChannel(targ.hdr,trigset);
    end
    targ.samplelimits = [1 targ.hdr.nSamples-targ.numsamples];
    targ.type = 'EDF';
else
    error('Unrecognized file format.');
end

% get target trigger edges and separations
targ.trigedges = formatTriggerEdges(targ.trigedges,targ.fs,trigset);
targ.trigseps = computeTriggerSeparations(targ.trigedges,targ.fs,trigset);

% initialize EDF sample limits
orig.fs = getFileSampleRate(originpath);
orig.numsamples = options.origintargetratio*round(options.searchlength*orig.fs);
if isNsxPath(originpath)
    warning('This has not been tested!');
    if isempty(orig.trigchan)
        orig.trigchan = getNsxTriggerChannel(originpath,trigset);
    end
    nevpath = regexprep(originpath,'.ns\d','.nev');
    nev = openNEV(nevpath,'nosave','nomat');
    orig.samplelimits = [1 nev.MetaTags.DataDuration-orig.numsamples];
    orig.type = 'NSx';
elseif isNevPath(originpath)
    % all NEV trigger events are preloaded
    nev = openNEV(originpath,'nosave','nomat');
    orig.alltrigedges = getNevTriggerEdges(nev);
    orig.samplelimits = [1 nev.MetaTags.DataDuration-orig.numsamples];
    orig.type = 'NEV';
elseif isEdfPath(originpath)
    orig.hdr = ft_read_header(originpath);
    if isempty(orig.trigchan)
        orig.trigchan = getEdfTriggerChannel(orig.hdr,trigset);
    end
    orig.samplelimits = [1 orig.hdr.nSamples-orig.numsamples];
    orig.type = 'EDF';
else
    error('Unrecognized file format.');
end

% if we know more, restrict our limit further
if ~isempty(options.lower), orig.samplelimits(1) = options.lower; end
if ~isempty(options.upper), orig.samplelimits(2) = options.upper; end

isfound = false;
numattempts = 0;
fprintf('\tSearching for %s triggers:\n',targ.type);
% do a binary search for the trigger block in the EDF file
while true
    % choose a new EDF search window
    if ~isfound
        midpoint = orig.samplelimits(1) + round(diff(orig.samplelimits)/2);
        orig.window = [midpoint, midpoint + orig.numsamples-1];
    end

    % load EDF triggers and get pulse edges/separations
    if isNsxPath(originpath)
        warning('This has not been tested!');
        nsx = openNSx(srchpath,'channels',orig.trigchan,'duration',...
            orig.window(1):orig.window(2));
        orig.data = nsx.Data;
        orig.trigedges = getAnalogTriggerEdges(orig.data);
    elseif isNevPath(originpath)
        % search the entire trigger pack
        orig.trigedges = orig.alltrigedges;
    elseif isEdfPath(originpath)
        orig.data = ft_read_data(originpath,'header',orig.hdr,'chanindx',...
            orig.trigchan,'begsample',orig.window(1),'endsample',orig.window(2));
        orig.trigedges = getAnalogTriggerEdges(orig.data);
    end
    orig.trigedges = formatTriggerEdges(orig.trigedges,orig.fs,trigset);
    orig.trigseps = computeTriggerSeparations(orig.trigedges,orig.fs,trigset);

    % try to match our target and origin trigger separations
    [hasoverlap,offsetidx] = matchTriggerSeparations(targ.trigseps,orig.trigseps);

    if hasoverlap
        % we've found at least part of the target trigger block
        % zero offset indicates that we've found the first and last trigger
        % otherwise shift the EDF sample window in the right direction and
        % repeat to capture all target triggers
        isfound = true;
        if ~offsetidx
            % get the edge index after getting the EDF separation group
            % use it to recover the timing
            [~,sepidx] = ismember(targ.trigseps(1,:),orig.trigseps,'rows');
            numpulsespergroup = size(targ.trigseps,2);
            pulseidx = 1+(sepidx-1)*numpulsespergroup;

            % sync the edges and leave the loop
            % note that for NEV searching the window isn't used
            if isNevPath(originpath)
                t0.origin = orig.trigedges(pulseidx,1);
            else
                t0.origin = orig.window(1) + orig.trigedges(pulseidx,1);
            end
            t0.target = targ.t0 + targ.trigedges(1,1);
            break
        else
            % shift the sample window to capture all triggers
            % try to match all the triggers in the next loop iteration
            if offsetidx > 0
                dt = 0.9*(targ.trigedges(offsetidx,1)-targ.trigedges(1,1))/targ.fs;
            else
                dt = -1.1*(targ.trigedges(-offsetidx,1)-targ.trigedges(1,1))/targ.fs;
            end
            orig.window = orig.window + round(dt*orig.fs);
            continue
        end
    else
        % compare the target and origin triggers to work out which direction to
        % shift the origin window
        shiftdir = getNextTriggerSearchDirection(targ.trigseps(1,:),orig.trigseps(1,:),trigset);
        if shiftdir > 0
            orig.samplelimits = [midpoint, orig.samplelimits(2)];
        else
            orig.samplelimits = [orig.samplelimits(1), midpoint];
        end

        fprintf('\t\t * %s window [%d-%d]\n',orig.type,...
            orig.samplelimits(1),orig.samplelimits(2));
        numattempts = numattempts + 1;
    end

    % quit if we're taking too long (maybe the target triggers aren't here)
    if numattempts >= options.maxattempts
        warning('Exceeded max attempts; could not find target trigger.');
        t0 = nan;
        return
    end
end

fprintf('\tFound first %s trigger (n=%d) in %s file (n=%d).\n',...
    targ.type,t0.target,orig.type,t0.origin);

if options.validate
    % for loading NEV events, allow a little buffer
    if isNevPath(targetpath) || isNevPath(originpath)
        nevpre = options.nevpresamples;
    end

    % load up matching data from target and EDF
    targ.window = [t0.target, t0.target+targ.numsamples-1];
    if isNsxPath(targetpath)
        nsx = openNSx(targetpath,'channels',targ.trigchan,'duration',...
            targ.window(1):targ.window(2));
        targ.data = nsx.Data;
    elseif isNevPath(targetpath)
        targ.nev = openNEV(targetpath,'nosave','nomat');
        targ.trigedges = getNevTriggerEdges(targ.nev,'window',targ.window-nevpre);
        targ.trigedges = formatTriggerEdges(targ.trigedges,targ.fs,trigset);
        [targ.time,targ.data] = makeWaveformFromEdges(targ.trigedges,'fs',targ.fs);
    elseif isEdfPath(targetpath)
        warning('This has not been tested!')
        targ.data = ft_read_data(targetpath,'header',targ.hdr,'chanindx',targ.trigchan,...
            'begsample',targ.window(1),'endsample',targ.window(2));
    end

    % convert time to seconds
    if ~isfield(targ,'time'), targ.time = (1:numel(targ.data))/targ.fs; end

    if ~isNevPath(targetpath)
        % clean up analog signal amplitude
        upthresh = options.maxamppc*max(targ.data);
        targ.data(targ.data <= upthresh) = 0;
        targ.data(targ.data > upthresh) = 1;
    end

    orig.numsamples = round(options.searchlength*orig.fs);
    orig.window = [t0.origin, t0.origin+2*orig.numsamples-1];
    if isNsxPath(originpath)
        warning('This has not been tested!');
        nsx = openNSx(originpath,'channels',orig.trigchan,'duration',...
            orig.window(1):orig.window(2));
        orig.data = nsx.Data;
    elseif isNevPath(originpath)
        orig.nev = openNEV(originpath,'nosave','nomat');
        orig.trigedges = getNevTriggerEdges(orig.nev,'window',orig.window-nevpre);
        orig.trigedges = formatTriggerEdges(orig.trigedges,orig.fs,trigset);
        [orig.time,orig.data] = makeWaveformFromEdges(orig.trigedges,'fs',targ.fs);
        orig.data = 0.5*orig.data;
    elseif isEdfPath(originpath)
        orig.data = ft_read_data(originpath,'header',orig.hdr,'chanindx',orig.trigchan,...
            'begsample',orig.window(1),'endsample',orig.window(2));
    end

    % convert time to seconds
    if ~isfield(orig,'time'), orig.time = (1:numel(orig.data))/orig.fs; end

    if ~isNevPath(originpath)
        % clean up analog signal amplitude
        upthresh = options.maxamppc*max(orig.data);
        orig.data(orig.data <= upthresh) = 0;
        orig.data(orig.data > upthresh) = 0.5;
    end

    % validation plot
    figure
    set(gcf,'color','w')
    hold on
    plot(targ.time,targ.data);
    plot(orig.time,orig.data);
    xlabel('seconds','fontsize',22)
    set(gca,'tickdir','out')
    set(gca,'ticklen',0.02*[1 1])
    legs = {sprintf('target %s',targ.type),sprintf('origin %s',orig.type)};
    legend(legs,'fontsize',16,'location','northeast')
end

end

% mini helper functions
function isnsx = isNsxPath(fpath)
isnsx = endsWith(fpath,{'.ns3','.ns5'});
end

function isnev = isNevPath(fpath)
isnev = endsWith(fpath,{'.nev'});
end

function isedf = isEdfPath(fpath)
isedf = endsWith(fpath,{'.edf'});
end
