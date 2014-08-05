function [meanR,meanV]=getRVmeansDDFTPlanar2D(ddft)
% [meanr,meanv]=getRVmeansDDFT(ddft,geom)
%   returns mean position and velocity for DDFT results at each time
%
% INPUTS:
%   ddft -- a structure containing
%           rhov        (2*N,nTimes) matrix with columns [rho; v] specified at 
%                           grid points r and each of the nTimes
%           r           (N,1) vector of grid points
%           w           (1,N) vector of corresponding weights
%           meanR       the mean position at each time step (optional)
%           meanP       the mean momentum at each time step (optional)
%           geom        geometry, should be 'planar' or 'spherical'
%
% OUTPUTS:
%   meanr       -- (nTimes,1) vector of mean positions
%   meanv       -- (nTimes,1) vector of mean velocities

% if the code doesn't compute the means, do it here

% assign relevant data
rho_t=ddft.dynamicsResult.rho_t;
y1=ddft.IDC.Pts.y1_kv;
y2=ddft.IDC.Pts.y2_kv;
Int=ddft.IDC.Int;

if(ndims(rho_t)==2)
    temp=zeros(size(rho_t,1),1,size(rho_t,2));
    temp(:,1,:)=rho_t;
    rho_t=temp;
end

% number of times in total
nTimes=size(rho_t,3);
% number of species
nSpecies=size(rho_t,2);
dim=2;

% mean vectors
meany1=zeros(nTimes,nSpecies);
meany2=zeros(nTimes,nSpecies);
%meanv=meanr;

y1(abs(y1)==Inf)=0;
y2(abs(y2)==Inf)=0;

% calculate mean at each time -- may be able to vectorize this
for iTime=1:nTimes
    % calculate rho and v from full rhov
    rho=rho_t(:,:,iTime);

    for iSpecies=1:nSpecies
        rhoS=rho(:,iSpecies);
        % calculate (normalised) means
        meany1(iTime,iSpecies)=( Int*(rhoS.*y1) ) / (Int*rhoS);
        meany2(iTime,iSpecies)=( Int*(rhoS.*y2) ) / (Int*rhoS);
    end
end

meanR=zeros(nTimes,nSpecies,dim);  

meanR(:,:,1)=meany1;
meanR(:,:,2)=meany2;

meanV=zeros(nTimes,nSpecies,dim);  

end