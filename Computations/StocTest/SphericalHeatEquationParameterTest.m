close all

Nvals = 10:10:200;
Lvals = 1:1:14;

maxRelErr = zeros(length(Nvals),length(Lvals));
t         = maxRelErr;
DErr      = maxRelErr;
D2Err     = maxRelErr;
DDErr     = maxRelErr;
IntErr    = maxRelErr;

for iN = 1:length(Nvals)
    for iL = 1:length(Lvals)
        [maxRelErr(iN,iL),t(iN,iL),DErr(iN,iL),D2Err(iN,iL),DDErr(iN,iL),IntErr(iN,iL)] ...
            = SphericalHeatEquation(Nvals(iN),Lvals(iL));
    end
end
