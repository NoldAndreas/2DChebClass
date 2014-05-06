function [R,P]=getRP(r,p,geom,dim)
%  [R,P]=getRP(r,p,geom,dim)
%   returns vectors of the appropriate one-dimensional positions for 
%   vectors x and momenta p in dimension dim and geometry geom.
%
% INPUTS:
%   x           -- (dim*nParticles,nSamples) matrix of positions 
%   p           -- (dim*nParticles,nSamples) matrix of momenta 
%   dim         -- dimension
%   geom        -- geometry, should be 'planar' or 'spherical'
%
% OUTPUTS:
%   R          -- (nParticles,dim,nSamples) matrix of dim-dimensional positions
%   P          -- (nParticles,dim,nSamples) matrix of dim-dimensional momenta


% get number of particles and number of samples
nParticles= size(r,1)/dim;

switch geom
    case 'planar'
        % the dim-th (i.e. last) row of mask contains the positions of the
        % dim-th coordinate of each particle
        mask=reshape(1:size(r,1),dim,nParticles);
        
        % get the dim-th coordinate of each particle for both position and
        % momentum
        R=r(mask(dim,:),:);
        P=p(mask(dim,:),:);

    case 'spherical'
        % reshape r and p so that each page contains the coordinates for
        % one sample with each column containing the coordinates for one
        % particle
        r=reshape(r,dim,nParticles,[]);
        p=reshape(p,dim,nParticles,[]);
        
        % calculate the norm of each column
        R=squeeze(sqrt(sum(r.^2,1)));
        % calculate r.p/|r|
        P=squeeze(sum(r.*p,1))./R;
        
    case 'planar2D'
        % get both coordinates - output is a matrix with two columns with
        % x and y for each particle in the rows
        
        % r, p are nParticles*dim x nSamples
        
        r=reshape(r,dim,nParticles,[]);
        p=reshape(p,dim,nParticles,[]);
        
        % r, p are now dim x nParticles x nSamples
        % with each particle's coordinates in a column and each sample on a
        % separate page
        
        % R, P should be nPartices x dim x nSamples so swap 1st and 2nd set
        % of indices
        R=permute(r,[2 1 3]);
        P=permute(p,[2 1 3]);

    case 'polar2D'

        % r, p are nParticles*dim x nSamples
        r=reshape(r,dim,nParticles,[]);
        p=reshape(p,dim,nParticles,[]);
         
        % r, p are now dim x nParticles x nSamples
        % with each particle's coordinates in a column and each sample on a
        % separate page
        
        % get polar coordinates for each particle and each sample
        [thetar,absr]=cart2pol(r(1,:,:),r(2,:,:));
        thetar=mod(thetar,2*pi);
        % these are size 1 x nParticles x nSamples
        
        % construct radial and tangential components of momentum
        pr     =  p(1,:,:).*cos(thetar) + p(2,:,:).*sin(thetar);
        ptheta = -p(1,:,:).*sin(thetar) + p(2,:,:).*cos(thetar);
        
        % reform into output matrices
        nSamples=size(r,3);
        R=zeros(nParticles,dim,nSamples);
        P=R;
        
        R(:,1,:)=permute(absr,[2 1 3]);
        R(:,2,:)=permute(thetar,[2 1 3]);
        
        P(:,1,:)=permute(pr,[2 1 3]);
        P(:,2,:)=permute(ptheta,[2 1 3]);

end

end %function


%         %check by plotting     
%         sample=1;
%         pr=pr(:,:,sample)';
%         ptheta=ptheta(:,:,sample)';
%         thetar=thetar(:,:,sample)';
%         r1=r(1,:,sample)';
%         r2=r(2,:,sample)';
%         p1=p(1,:,sample)';
%         p2=p(2,:,sample)';
%         
%         rhat=[cos(thetar) sin(thetar)];
%         
%         prPlot=[pr pr].*rhat + [r1 r2];
%         
%         thetahat = [-sin(thetar) cos(thetar)];
%         
%         pthetaPlot= [ptheta ptheta].*thetahat + prPlot;
%         
%         for iParticle=1:nParticles
%             plot([0 r1(iParticle)],[0 r2(iParticle)],'k');
% 
%             hold on;
%             plot([r1(iParticle) r1(iParticle)+p1(iParticle)],[r2(iParticle) r2(iParticle)+p2(iParticle)],'r');
%             
%             plot([r1(iParticle) prPlot(iParticle,1)],[r2(iParticle) prPlot(iParticle,2)],'--b');
%             
%             plot([prPlot(iParticle,1) pthetaPlot(iParticle,1)],[prPlot(iParticle,2) pthetaPlot(iParticle,2)],'--g');
%             
%         end
%         axis 'equal';