function [shiftdir] = getNextTriggerSearchDirection(targetseps, sourceseps, trigset)
% GETNEXTTRIGGERSEARCHDIRECTION returns the direction (+1 to shift the
%	target window up, -1 to shift it down) of the shift in the next trigger
%	search window. A shift of 0 means we are probably confused...
%
% Alex Hadjinicolaou <a.e.hadjinicolaou@gmail.com>

if strcmp(trigset,'presentation')
    % first and second separations encode hours and minutes
    shiftdir = sign(targetseps(1)-sourceseps(1));
    if ~shiftdir
        shiftdir = sign(targetseps(2)-sourceseps(2));
    end
else
    error('Unrecognized trigset (%s).',trigset);
end
