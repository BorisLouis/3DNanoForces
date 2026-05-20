% --- USER SETTINGS ---
tic;
InputDirectory = 'F:\WoutVT\20251023_200nm_PS_2\20251023_200nm_PS_2\1MHz_3V';
load([InputDirectory '\' 'trackResults.mat']);

outFile = '200nm_1MHz_3V.h5';
expTime = 0.008;
% --- DELETE IF EXISTS (avoid append confusion) ---
if exist(outFile,'file'); delete(outFile); end
% --- FIRST PASS: count total rows (points) ---
nTraces = size(trackRes.traces, 1);
totalN = 0;
for id = 1:nTraces
    currTrace = trackRes.traces{id,1};
    if isempty(currTrace); continue; end
    totalN = totalN + numel(currTrace.t);
end

% --- CREATE DATASET ---
% Columns: [id, t, x, y, z]
h5create(outFile, '/tracks', [totalN 5], 'Datatype', 'double', 'ChunkSize', [min(totalN,1e5) 5]);

% Optional: store column names + units as attributes
h5writeatt(outFile, '/tracks', 'columns', 'id,t,x,y,z');
h5writeatt(outFile, '/tracks', 'units',   'id:arb; t:s; x,y:um (or your unit); z:yourunit');

% --- SECOND PASS: fill data ---
k = 1;
for id = 1:nTraces
    currTrace = trackRes.traces{id,1};
    if isempty(currTrace); continue; end

    x = currTrace.col ;
    y = currTrace.row ;
    z = currTrace.z;
    t   = (currTrace.t-currTrace.t(1))*expTime;
    t   = t/(max(t));

    n = numel(t);
    block = [ ...
        repmat(double(id), n, 1), ...
        double(t(:)), ...
        double(x(:)), ...
        double(y(:)), ...
        double(z(:)) ...
    ];

    h5write(outFile, '/tracks', block, [k 1], [n 5]);
    k = k + n;
end

fprintf('Wrote %d points to %s\n', totalN, outFile);

toc;
