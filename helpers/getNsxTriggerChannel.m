function [idx] = getNsxTriggerChannel(nsxpath, trigset)
% GETNSXTRIGGERCHANNEL Gets the trigger channel for the input NSx file.
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

switch trigset
    case 'presentation'
        triglabel = 'image';
    otherwise
        error('Unrecognized trigset (%s).',trigset);
end

% read a sample from the NSx file and search for the trigger
nsx = openNSx(nsxpath,'t:1:2');
idx = find(strcmpi(deblank({nsx.ElectrodesInfo.Label}),triglabel));

if isempty(idx)
    warning('Channel %s not found in NSx file.',triglabel);
elseif numel(idx) > 1
    warning('Channel %s is not unique in NSx file.',triglabel);
end
