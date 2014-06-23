function [nR,meanP,xR]=RPhistShift(R,P,nBins,Rmin,Rmax,shift)
%[nR,nP,xR]=RPhist(R,meanP,nBins)
%   determines histogram data for positions and mean momentum in each of
%   the histogram bins
%
% INPUTS: 
%  R            -- (nSamples*nParticles,1) vector of 1D positions
%  P            -- (nSamples*nParticles,1) vector of 1D momenta
%  nBins        -- number of bins for histogram
%
% OUTPUTS:
%  nR            -- (nBins,1) vector of the number of particles in each bin
%  meanP         -- (nBins,1) vector of the mean momentum of each bin
%  xR            -- (nBins,1) vector of bin centers

if(nargin<4)
    % find bounds on R data, and add a bit of a buffer to ensure we capture
    % everything
    Rmax=max(R)+10*eps;
    Rmin=min(R)-10*eps;
    shift=0;
end

% calculate (uniform) bin width
binWidth=(Rmax-Rmin)/nBins;

% calculate ends of bins
binEnds=Rmin:binWidth:Rmax; % has length nBins+1

% add buffer bins at the end
binEnds=[Rmin-binWidth, binEnds, Rmax+binWidth];
nBins=nBins+2;

binEnds=binEnds+shift;

% set up output vectors
nR=zeros(nBins,1);
meanP=nR;

for iBin=1:nBins  % loop through bins
    % find which R are in the bin
    binMask=(binEnds(iBin)<=R & R < binEnds(iBin+1));
    % count how many
    nR(iBin)=sum(binMask);
    % find mean momentum
    if(nR(iBin)>0)
        meanP(iBin)=sum(P(binMask))/nR(iBin);
    end
end

% find centres of bins
xR=binEnds+0.5*binWidth;
xR=xR(1:nBins);
xR=xR(:);

end

