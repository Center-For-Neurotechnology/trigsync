function [hasoverlap, offsetidx] = matchTriggerSeparations(targetseps, originseps)
% MATCHTRIGGERSEPARATIONS indicates whether the target trigger
%   pulse separations can be found within the origin separations structure.
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

% overlap: either the first pulse separations and/or the last ones can
% be found within the origin pulse separations
[firstfound,firstidx] = ismember(targetseps(1,:),originseps,'rows');
if firstfound, originseps(firstidx,:) = nan; end
[lastfound,lastidx] = ismember(targetseps(end,:),originseps,'rows');

hasoverlap = firstfound || lastfound;

offsetidx = nan;
if firstfound && lastfound
    offsetidx = 0;
elseif firstfound
    offsetidx = firstidx;
elseif lastfound
    offsetidx = -lastidx;
end
