function [meanr,meanv]=getRVmeansDDFT(ddft,geom)
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
if(~isfield(ddft,'meanR'))
    % NOT DONE FOR N COMPONENTS
    fprintf(1,'NOT DONE FOR N COMPONENTS');
    
    % assign relevant data
    rhov=ddft.rhov;
    r=ddft.r;
    w=ddft.w;

    % number of times in total
    %nTimes=size(rhov,1);
    nTimes=size(rhov,2);

    % mean vectors
    meanr=zeros(nTimes,1);
    meanv=meanr;

    % number of points in r
    N=size(rhov,1)/2;

    % and find appropriate integration weights
    switch geom    
        case 'spherical'
            weight=4*pi*r.^2.*w';
        case 'planar'
            weight=w';
    end

    % calculate mean at each time -- may be able to vectorize this
    for iTime=1:nTimes
        % calculate rho and v from full rhov
        rho=rhov(1:N,iTime);
        v=rhov(N+1:end,iTime);
        
        % calculate (normalised) means
        meanr(iTime,1)=sum(rho.*r.*weight)/sum(rho.*weight);
        meanv(iTime,1)=sum(rho.*v.*weight)/sum(rho.*weight);
    end

else
    % computed in the DDFT code so just read it it
    meanr=(ddft.meanR)';
    meanv=(ddft.meanP)';
end
    
end