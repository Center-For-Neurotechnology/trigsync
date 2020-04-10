function [idx] = getEdfTriggerChannel(hdr, trigset)
% GETEDFTRIGGERCHANNELS returns the index of a trigger channel in an EDF.
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

% allow hdr to be a path to the header file (instead of the actual header)
if ~isstruct(hdr)
    if endsWith(hdr,'.edf')
        hdr = ft_read_header(hdr);
    else
        error('HDR input must be either an EDF header or a path to the header file.');
    end
end

switch trigset
    case 'presentation'
        triglabel = 'trig';
    otherwise
        error('Unrecognized trigset (%s).',trigset);
end

idx = find(strcmpi(deblank(hdr.label),triglabel));

if isempty(idx)
    warning('Channel %s not found in EDF file.',triglabel);
elseif numel(idx) > 1
    warning('Channel %s is not unique in EDF file.',triglabel);
end
