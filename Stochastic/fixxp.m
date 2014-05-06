function [x,p,percentage]=fixxp(x,p,xLower,xUpper)
%[xCut,pCut,percentage]=stripxp(x,p,xLower,xUpper)
%  replace calculations where the time step was too large or the
%  potential too steep, causing particles to fly off to infinity
% 
% INPUTS
%  x           --  positions (nRuns,nParticles,nTimes)
%  p           --  momenta   (nRuns,nParticles,nTimes)
%  xLower      --  lower bound to cut off beyond
%  xUpper      --  upper bound to cut off beyond
%
% OUTPUTS
%  x           --  amended positions (?,nParticles,nTimes)
%  p           --  amended momenta   (?,nParticles,nTimes)
%  percentage  --  number of runs elminated

% get number of runs
nRuns=size(x,1);

% default mask is to replace nothing
mask=false(nRuns,1);

% number of runs we've cut out
nCut=0;

MAXPos = 0;
MINPos = 0;

% check each run
for iRun=1:nRuns
    % determine if any particles are outside the cutoff
    minPos=min( min( x(iRun,:,:) ) );
    maxPos=max( max( x(iRun,:,:) ) );
    
    MAXPos = max(MAXPos,maxPos);
    MINPos = min(MINPos,minPos);
    
    outOfBox=(maxPos > xUpper || isnan(maxPos) || minPos < xLower || isnan(minPos));
    if( outOfBox )
        % if so then eliminate this run
        mask(iRun)=true;
        % add one to the number we've eliminated
        nCut=nCut+1;
    end
end

MAXPos
MINPos

% replace by first correct value
replacePos=find(~mask,1,'first');
xReplace=x(replacePos,:,:);
xReplace=repmat(xReplace,nCut,1);
pReplace=p(replacePos,:,:);
pReplace=repmat(pReplace,nCut,1);

% set output
x(mask,:,:)=xReplace;
p(mask,:,:)=pReplace;



% and percentage we've cut out
percentage=nCut/nRuns;