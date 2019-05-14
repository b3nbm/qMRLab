function [conc, signal0] = ConcFromSignal_SPGR(data, flipAngle, TR, T1Map, relaxivity1, b1Map, roi) %Compute CA concentration of dynamic contrast enhanced SPGR data 
%
% function conc = ConcFromSignal_SPGR(data, flipAngle, TR, T1Map, roi, relaxivity1, b1Map)
% -----------------------------------------------------------
% INPUTS:
%   data: width x length x slices x flipAngle matrix
%   flipAngle: flip angle (in degrees) 
%   TR: Repetition time in seconds
%   T1Map: width x length x slices matrix containing T1 (s)
%   roi: width x length x slices binary mask 
%   reaxivity1: CA R1 relaxivity (mM/ms)
%   b1Map: width x length x slices matrix containing relative flip angle
%         (i.e. if nominal alpha is 60 and measured alpha is 61, then b1Map = 61/60
%   verbose: logical - for debugging


if ndims(data)<3, data = permute(data(:),[2 3 4 1]); end
[nX, nY, nZ, nTime] = size(data);
nVox = nX*nY*nZ;

if ~exist('b1Map', 'var') || isempty(b1Map)
    b1Map = ones(nX, nY, nZ);
end

if ~exist('roi', 'var') || isempty(roi)
    roi = true(nX, nY, nZ);
end

if ~exist('verbose', 'var')
    verbose = 0;
end
 
if length(b1Map) ~= length(data(:,:,:,1)), error('B1 size is different from data size'); end
if length(T1Map) ~= length(data(:,:,:,1)), error('T1 size is different from data size'); end
if ~islogical(roi), roi = logical(roi); end

% Reshape data into 2D array so that future steps are simpler
data = reshape(data,[nVox nTime])'; % Transpose because MATLAB is column-major
data = data(:, roi(:));
nVox = sum(roi(:));

% Reshape T1Map and B11Map into 1D array so that future steps are simpler
R10 = 1./T1Map(roi(:))';
b1Map = b1Map(roi(:))';

% Large matrix with redundant values will make computations faster,
% but will also consume more memory
alpha = deg2rad(flipAngle);
if isrow(alpha), alpha = alpha'; end
alpha = repmat(alpha,[1 nVox]) .* b1Map;

% Compute Dynamic Concentration Data
warning('Initial signal (S0) is set to first frame.');
signal0 = data(1,:);
dr1 = signal2relax(data,R10,signal0,TR,alpha);

[conc] = zeros(nX,nY,nZ,nTime);
for iTime = 1:nTime
    currentConc = zeros(nX,nY,nZ);
    currentConc(roi(:)) = dr1(iTime,:)/relaxivity1;
    conc(:,:,:,iTime) = currentConc;
end

end % END OF DCE_OnSPGR

function dr1 = signal2relax(signal,R10,signal0,TR,alpha)
% COMPUTE DATA CONC
f10 = (1-cos(alpha).*exp(-TR.*R10)) ./ (1-exp(-TR.*R10));
[nTime,nVox] = size(signal);
dr1 = NaN(nTime,nVox);
for iT=1:nTime
    currentSignalRatio = signal(iT,:)./signal0;
    currentDr1 = log((f10-cos(alpha).*currentSignalRatio) ./ ...
        (f10 - currentSignalRatio) )/TR - R10;
    dr1(iT,:) = currentDr1;
end

end
