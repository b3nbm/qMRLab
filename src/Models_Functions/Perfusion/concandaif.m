function [conc, aif] = concandaif(data, flipAngle, TR, T1Map, aifMask, relaxivity1, b1Map) %Perform a simple non linear-least squares data fit on dynamic contrast enhanced SPGR data 
%
% function [KTrans, Ve, Vb, AIF] = DCE_OnSPGR(data, time, flipAngle, TR, T1Map, aifMask, hematocrit, relaxivity1, b1Map, roi, verbose)
% -----------------------------------------------------------
% INPUTS:
%   data: width x length x slices x flipAngle matrix
%   time: frame time (in seconds) corresponding to last dimension of 'data'
%   flipAngle: flip angle (in degrees) 
%   TR: Repetition time in seconds
%   b1Map: width x length x slices matrix containing relative flip angle
%         (i.e. if nominal alpha is 60 and measured alpha is 61, then b1Map = 61/60
%   roi: width x length x slices binary mask 
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
if ~islogical(aifMask), aifMask = logical(aifMask); end
if ~islogical(roi), roi = logical(roi); end

% Reshape data into 2D array so that future steps are simpler
data = reshape(data,[nVox nTime])'; % Transpose because MATLAB is column-major
data = data(:, roi(:));
nVox = sum(roi(:));

% Reshape T1Map into 1D array so that future steps are simpler
R10 = 1./T1Map(roi(:))';

% Large matrix with redundant values will make computations faster,
% but will also consume more memory
alpha = deg2rad(flipAngle);
if isrow(alpha), alpha = alpha'; end
alpha = repmat(alpha,[1 nVox]) .* b1Map(roi(:))';

% Compute Dynamic Concentration Data
signal0 = data(1,:);
conc = signal2conc(data,relaxivity1,R10,signal0,TR,alpha);

% Compute AIF from mask
if ~exist('aifMask', 'var') || isempty(aifMask)
    aif = [];
else
    aifAllConc = conc(:, aifMask(:));
    aif = mean(aifAllConc,2);
end

end % END OF DCE_OnSPGR

function conc = signal2conc(signal,r1,R10,signal0,TR,alpha)
% COMPUTE DATA CONC
f10 = (1-cos(alpha).*exp(-TR.*R10)) ./ (1-exp(-TR.*R10));
[nT,nVox] = size(signal);
dr1 = NaN(nT,nVox);
for iT=1:nT
    currentSignalRatio = signal(iT,:)./signal0;
    currentDr1 = log((f10-cos(alpha).*currentSignalRatio) ./ ...
        (f10 - currentSignalRatio) )/TR - R10;
    dr1(iT,:) = currentDr1;
end
conc = dr1/r1;
end

function conc = pkm(x,time,aif)
% Extended Tofts model
% --- Convolution method adapted for variable time vectors
% Evaluate IRF model for all ca time repeated for each ct time values
KTrans = x(1);
Ve = x(2);
Vb = x(3);
kep = KTrans/Ve;

dt = time(end)-time(end-1);
dtVar = [(time(2:end) - time(1:(end-1))); dt];
if any(dtVar - dt > 0.001)
    warning('Convolution fails with variable time sampling. Although I don''t know exactly why... should be investigated.')
end
[tauM, timeM] = meshgrid(time, time);
timeConvDiffM = timeM - tauM;

biexpIrf = exp(-kep*timeConvDiffM); % (-)
biexpIrf(timeConvDiffM<0) = 0;

conc = biexpIrf*(aif.*dtVar)+Vb*aif;

conc = KTrans*conc;
end
