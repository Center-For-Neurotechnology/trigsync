function [fs] = getFileSampleRate(fpath)
% GETFILESAMPLERATE Gets the sample rate by working out what kind of data
%   the filepath represents.
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

if isNsxPath(fpath)
    % Blackrock NSx
    nsx = openNSx(fpath,'duration',1:2);
    fs = nsx.MetaTags.SamplingFreq;
elseif isNevPath(fpath)
    % Blackrock NEV
    nev = openNEV(fpath,'nosave','nomat');
    fs = double(nev.MetaTags.SampleRes);
elseif isEdfPath(fpath)
    % clinical EDF
    ehdr = ft_read_header(fpath);
    fs = ehdr.Fs;
else
    error('Unrecognized file format.');
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