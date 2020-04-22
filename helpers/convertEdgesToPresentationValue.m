function hms = convertEdgesToPresentationValue(edges,fs)

% complain if not enough elements
if numel(edges) < 7, error('Not enough edges.'); end

% flatten into array
if size(edges,1) > 1
    edges = edges';
    edges = edges(:);
end

% only use the first 7 edges
if numel(edges) > 7, edges = edges(1:7); end

% convert to time separations
hms = [...
    edges(3)-edges(2),...
    edges(5)-edges(4),...
    edges(7)-edges(6)]/fs;

% round to units of 50ms
roundunit = 5e-2;
hms = round(hms/roundunit)*roundunit;
