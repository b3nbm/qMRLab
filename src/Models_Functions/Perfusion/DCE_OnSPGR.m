function [KTrans, Ve, Vb, aif] = DCE_OnSPGR(conc, time, aif, st, lb, ub, fx, brainDensity, verbose)
%Perform a simple non linear-least squares data fit on dynamic contrast enhanced SPGR data 
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


if ndims(conc)<3, conc = permute(conc(:),[2 3 4 1]); end
[nX, nY, nZ, nTime] = size(conc);
nVox = nX*nY*nZ;

if ~exist('roi', 'var') || isempty(roi)
    roi = true(nX, nY, nZ);
end

if ~exist('verbose', 'var')
    verbose = 0;
end
 
if ~islogical(roi), roi = logical(roi); end

% Reshape data into 2D array so that future steps are simpler
conc = reshape(conc,[nVox nTime])'; % Transpose because MATLAB is column-major
conc = conc(:, roi(:));
nVox = sum(roi(:));

% Pharmacokinetic model
estKTrans = NaN(1,nVox);
estVe = NaN(1,nVox);
estVb = NaN(1,nVox);
lsqeq = @(x,xdata) PKM_eTofts(x,xdata,aif);
for iVox=1:nVox
    [xopt, resnorm] = lsqcurvefit(lsqeq,...
        st(~fx), time, conc(:,iVox), lb(~fx), ub(~fx));
    estKTrans(iVox) = xopt(1);
    estVe(iVox) = xopt(2);
    estVb(iVox) = xopt(3);
end
% Assign arbitrary KTrans, Ve and Vb value if fitted value is unphysical. 
% Might be better to set this to NaN so users can identify voxels where fit fails
failedFit = isnan(estKTrans) | isnan(estVb)  | isnan(estVe) ;
failedFitValue = NaN;
estKTrans(failedFit) = failedFitValue;
estVe(failedFit) = failedFitValue;
estVb(failedFit) = failedFitValue;

% Assign estimated KTrans, Ve and Vb values to correct voxel 
[KTrans, Ve, Vb] = deal(zeros(nX,nY,nZ));
KTrans(roi(:)) = estKTrans;
Ve(roi(:)) = estVe;
Vb(roi(:)) = estVb;
end % END OF DCE_OnSPGR

