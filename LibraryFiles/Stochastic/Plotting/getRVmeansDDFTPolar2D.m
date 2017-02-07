function [meanR,meanV]=getRVmeansDDFTPolar2D(ddft)
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
rho_t=ddft.rho_t;
y1=ddft.Pts.y1_kv;
y2=ddft.Pts.y2_kv;
Int=ddft.Int;

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
meany1C=zeros(nTimes,nSpecies);
meany2C=zeros(nTimes,nSpecies);
%meanv=meanr;

y1(y1==Inf)=0;
y2(y2==Inf)=0;

% calculate mean at each time -- may be able to vectorize this
for iTime=1:nTimes
    % calculate rho and v from full rhov
    rho=rho_t(:,:,iTime);

    for iSpecies=1:nSpecies
        rhoS=rho(:,iSpecies);
        % calculate (normalised) means
        meany1C(iTime,iSpecies)=( Int*(rhoS.*y1.*cos(y2)) ) / (Int*rhoS);
        meany2C(iTime,iSpecies)=( Int*(rhoS.*y1.*sin(y2)) ) / (Int*rhoS);
        
%         [~,maxPos]=max(rhoS);
%         
%         meany1C(iTime,iSpecies)=y1(maxPos)*cos(y2(maxPos));
%         meany2C(iTime,iSpecies)=y1(maxPos)*sin(y2(maxPos));
    end
    
end


% mean vectors
meany1=zeros(nTimes,nSpecies);
meany2=zeros(nTimes,nSpecies);

for iSpecies=1:nSpecies

    [meany2(:,iSpecies),meany1(:,iSpecies)]=cart2pol(meany1C(:,iSpecies),meany2C(:,iSpecies));
    meany2(:,iSpecies)=mod(meany2(:,iSpecies),2*pi);

end

meanR=zeros(nTimes,nSpecies,dim);  

meanR(:,:,1)=meany1;
meanR(:,:,2)=meany2;

% meanR(:,:,1)=meany1C;
% meanR(:,:,2)=meany2C;
% figure
% plot(meany1C(:,1),meany2C(:,1))
% pause

meanV=zeros(nTimes,nSpecies,dim);  

end