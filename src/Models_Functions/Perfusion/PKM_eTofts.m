function conc = PKM_eTofts(x,time,aif)
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
