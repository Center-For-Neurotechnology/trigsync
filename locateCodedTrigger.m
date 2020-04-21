function [t0] = locateCodedTrigger(targetpath, originpath, trigset, varargin)
% LOCATECODEDTRIGGER returns the sample of the first trigger in a TARGET
%   ephys file, and the corresponding sample of that first trigger within
%   an ORIGIN ephys file. These triggers are coded, meaning that their
%   value is related to their location in the recording.
%
%   TARGETPATH and ORIGINPATH can point to NSx, NEV, or EDF data.
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
    'targetupthreshpc',[],...
    'origintrigchan',[],...
    'originupthreshpc',[],...
    'origintargetratio',2,...
    'getorigintrigbounds',false,...
    'figposition',[402 457 1100 300],...
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

t0 = [];
% override default trigger channels if user-specified
targ.trigchan = options.targettrigchan;
orig.trigchan = options.origintrigchan;

% initialize key structures
targ.name = 'target'; targ.path = targetpath;
orig.name = 'origin'; orig.path = originpath;

% allow for direct NEV/NSx structure input
for tt = [targ,orig]
    if isstruct(tt.path) && isfield(tt.path,'MetaTags')
        ext = tt.path.MetaTags.FileExt;
        tpath = fullfile(tt.path.MetaTags.FilePath,...
            sprintf('%s%s',tt.path.MetaTags.Filename,ext));
        if strcmp(tt.name,'target')
            if isNsxPath(ext), tt.nsx = targetpath; else, tt.nev = targetpath; end
        else
            if isNsxPath(ext), tt.nsx = originpath; else, tt.nev = originpath; end
        end
        tt.path = tpath;
    end
    
    checkFileCompatibility(tt.path,tt.name);
    tt.type = getFileDescription(tt.path);
    tt.fs = getFileSampleRate(tt.path);
    tt.windowlimit = false(1,2);
    if isempty(tt.trigchan)
        tt.trigchan = getTriggerChannel(tt.path,trigset);
    end
    
    % store the new trigger data
    if strcmp(tt.name,'target'), targ = tt; else, orig = tt; end
end

% start sample for target trigger window
targ.t0 = 1;
if ~isempty(options.targetoffsettime)
    targ.t0 = 1+floor(options.targetoffsettime*targ.fs);
end

% prepare target trigger window
targ.numsamples = round(options.searchlength*targ.fs);
targ.window = [targ.t0, targ.t0+targ.numsamples-1];

% get target trigger edges and separations
targ = getTriggerEdgeData(targ,trigset,options.targetupthreshpc);
targ.trigseps = computeTriggerSeparations(targ.trigedges,targ.fs,trigset);

% preload all trigger edges for NEV origin (not necessary?)
% preload header structure for EDF origin
if isNevPath(orig.path) && false
    nev = openNEV(orig.path,'nosave','nomat');
    orig.alltrigedges = getNevTriggerEdges(nev);
elseif isEdfPath(orig.path)
    orig.hdr = ft_read_header(orig.path);
end

% initialize origin sample limits
orig.numsamples = options.origintargetratio*round(options.searchlength*orig.fs);
orig.samplelimits = getFileSampleLimits(orig);

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
        orig.window = [midpoint, midpoint+orig.numsamples-1];
    end

    % load EDF triggers and get pulse edges/separations
    orig = getTriggerEdgeData(orig,trigset,options.originupthreshpc);    
    orig.trigseps = computeTriggerSeparations(orig.trigedges,orig.fs,trigset);

    % check to see if the loaded triggers differ from the last ones
    % if not, the target triggers are likely not present in the origin
    if isfield(orig,'prevseps')
        if isAlmostEqual(orig.trigseps(1,:),orig.prevseps(1,:)) && ...
                isAlmostEqual(orig.trigseps(end,:),orig.prevseps(2,:))
            warning('Target triggers are (likely) not in the origin recording.');
            return
        end
    end
    
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
            [~,sep1idx] = ismember(targ.trigseps(1,:),orig.trigseps,'rows');
            [~,sep2idx] = ismember(targ.trigseps(2,:),orig.trigseps,'rows');
            
            % adjust separation index to compensate for repeats
            if sep2idx-sep1idx > 1
                sep1idx = sep1idx + 1;
            end
                
            numpulsespergroup = size(targ.trigseps,2);
            pulseidx = 1+(sep1idx-1)*numpulsespergroup;

            % sync the edges and leave the loop
            % note that for NEV searching the window isn't used
            if isNevPath(orig.path)
                t0.origin = orig.trigedges(pulseidx,1);
            else
                t0.origin = orig.window(1) + orig.trigedges(pulseidx,1);
            end
            t0.target = targ.t0 + targ.trigedges(1,1);
            break
        else
            % if we're at a window limit and still haven't found the first
            % trigger, report the first overlapping trigger
            if any(orig.windowlimit)
                % TODO: sort this out
                warning('not yet handled -- need to find a good test case')
                break
            end
            
            % shift the sample window to capture all triggers
            % try to match all the triggers in the next loop iteration
            if offsetidx > 0
                dt = 0.9*(targ.trigedges(offsetidx,1)-targ.trigedges(1,1))/targ.fs;
            else
                dt = -1.1*(targ.trigedges(-offsetidx,1)-targ.trigedges(1,1))/targ.fs;
            end
            nextwindow = orig.window + round(dt*orig.fs);
            
            % make sure the next origin window is bounded
            % flag if we've reached the origin window limits
            orig.windowlimit(1) = nextwindow(1) <= 1;
            orig.windowlimit(2) = nextwindow(2) >= orig.samplelimits(2);
            if orig.windowlimit(1), nextwindow(1) = 1; end
            if orig.windowlimit(2), nextwindow(2) = orig.samplelimits(2); end
            orig.window = nextwindow;
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

        % keep track of the first/last trigger separations
        orig.prevseps = [orig.trigseps(1,:); orig.trigseps(end,:)];
        
        fprintf('\t\t * %s window [%d-%d]\n',orig.type,...
            orig.samplelimits(1),orig.samplelimits(2));
        numattempts = numattempts + 1;
    end

    % quit if we're taking too long
    if numattempts >= options.maxattempts
        warning('Exceeded max attempts; could not find target trigger.');
        return
    end
end

fprintf('\tFound first %s trigger (n=%d) in %s file (n=%d).\n',...
    targ.type,t0.target,orig.type,t0.origin);

if options.validate
    % collect the data for both target and origin, aligned to t0
    toffset = 1;
    for tc = {targ, orig}
        tt = tc{1};
        tt.numsamples = round(options.searchlength*tt.fs);
        tt.window = t0.(tt.name) + [0, tt.numsamples-1] - round(toffset*tt.fs);
        tt = getTriggerWaveData(tt,trigset);
        if strcmp(tt.name,'target'), targ = tt; else, orig = tt; end
    end

    tidx = t0.target-targ.window(1);
    oidx = t0.origin-orig.window(1);
    colors = lines(2);
    
    % validation plot
    figure
    set(gcf,'color','w')
    set(gcf,'position',options.figposition);
    hold on
    set(gca,'fontsize',16)
    plot(targ.time,targ.data,'color',colors(1,:));
    plot(orig.time,orig.data,'color',colors(2,:));

    % marker for first trigger pair
    if ~isNevPath(targ.path)
        plot(targ.time(tidx),targ.data(tidx),'o','color',colors(1,:),'markerfacecolor','w');
    end
    if ~isNevPath(targ.path)
        plot(orig.time(oidx),orig.data(oidx),'o','color',colors(2,:),'markerfacecolor','w');
    end
    xlabel('seconds','fontsize',22)
    set(gca,'tickdir','out')
    set(gca,'ticklen',0.012*[1 1])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    zoom xon
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

function checkFileCompatibility(fpath,ttype)
if ~(isNsxPath(fpath) || isNevPath(fpath) || isEdfPath(fpath))
    error('Incompatible %s file type.',ttype);
end
end

function rdesc = getFileDescription(fpath)
rdesc = 'unknown';
if isNsxPath(fpath)
    rdesc = 'NSx';
elseif isNevPath(fpath)
    rdesc = 'NEV';
elseif isEdfPath(fpath)
    rdesc = 'EDF';
end
end

function chan = getTriggerChannel(fpath,trigset)
chan = [];
if isNsxPath(fpath)
    chan = getNsxTriggerChannel(fpath,trigset);
elseif isNevPath(fpath)
    % should something be done here in the future?
elseif isEdfPath(fpath)
    chan = getEdfTriggerChannel(fpath,trigset);
end
end

function s = getTriggerEdgeData(s,trigset,upthreshpc)
if isNsxPath(s.path)
    nsx = openNSx(s.path,'channels',s.trigchan,'duration',...
        s.window(1):s.window(2));
    s.data = nsx.Data;
    s.trigedges = getAnalogTriggerEdges(s.data,'trigset',trigset,'fs',s.fs);
elseif isNevPath(s.path)
    % load the header if not present
    if ~isfield(s,'nev'), s.nev = openNEV(s.path,'nosave','nomat'); end
    s.trigedges = getNevTriggerEdges(s.nev,'window',s.window,'upthreshpc',upthreshpc);
elseif isEdfPath(s.path)
    % load the header if not present
    if ~isfield(s,'hdr'), s.hdr = ft_read_header(s.path); end
    s.data = ft_read_data(s.path,'header',s.hdr,'chanindx',...
        s.trigchan,'begsample',s.window(1),'endsample',s.window(2));
    s.trigedges = getAnalogTriggerEdges(s.data,'trigset',trigset,'fs',s.fs);
end
s.trigedges = formatTriggerEdges(s.trigedges,s.fs,trigset);
end

function s = getTriggerWaveData(s,trigset)
    % reconstructed waveform amplitude
    upval = 1; if strcmp(s.name,'origin'), upval = 0.5; end
    
    % load up matching data from target and EDF
    if isNsxPath(s.path)
        nsx = openNSx(s.path,'channels',s.trigchan,'duration',...
            s.window(1):s.window(2));
        s.data = nsx.Data;
    elseif isNevPath(s.path)
        % allow a little sample buffer for NEV events
        nevpre = 1e3;
        if ~isfield(s,'nev'), s.nev = openNEV(s.path,'nosave','nomat'); end
        s.trigedges = getNevTriggerEdges(s.nev,'window',s.window-nevpre);
        s.trigedges = formatTriggerEdges(s.trigedges,s.fs,trigset);
        [s.time,s.data] = makeWaveformFromEdges(s.trigedges,'fs',s.fs,'upval',upval);
    elseif isEdfPath(s.path)
        s.data = ft_read_data(s.path,'header',s.hdr,'chanindx',s.trigchan,...
            'begsample',s.window(1),'endsample',s.window(2));
    end

    % convert time to seconds
    if ~isfield(s,'time'), s.time = (1:numel(s.data))/s.fs; end

    if ~isNevPath(s.path)
        % clean up analog signal amplitude
        upthresh = 0.8*max(s.data);
        s.data(s.data < upthresh) = 0;
        s.data(s.data >= upthresh) = upval;
    end
end

function slims = getFileSampleLimits(s)
if isNsxPath(s.path)
    nevpath = regexprep(s.path,'.ns\d','.nev');
    nev = openNEV(nevpath,'nosave','nomat');
    slims = [1 nev.MetaTags.DataDuration-s.numsamples];
elseif isNevPath(s.path)
    nev = openNEV(s.path,'nosave','nomat');
    slims = [1 nev.MetaTags.DataDuration-s.numsamples];
elseif isEdfPath(s.path)
    slims = [1 s.hdr.nSamples-s.numsamples];
end
end
