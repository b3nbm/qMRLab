function aif = AIFFromMask(conc,aifMask)

if ndims(conc)<3, conc = permute(conc(:),[2 3 4 1]); end
[nX, nY, nZ, nTime] = size(conc);
nVox = nX*nY*nZ;

if ~islogical(aifMask), aifMask = logical(aifMask); end

conc = reshape(conc,[nVox nTime])'; % Transpose because MATLAB is column-major

% Compute AIF from mask
if ~exist('aifMask', 'var') || isempty(aifMask)
    aif = [];
else
    aifAllConc = conc(:, aifMask(:));
    aif = mean(aifAllConc,2);
end
