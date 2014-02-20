function [V1,DV1,VBack,DVBack,VAdd]=polar2Dtest(r,theta,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

switch optsPhys.type
    
    case 'stoc'
        
        % ensure column vector, e.g. for slice sampling which uses row
        % vectors
        x=r(:);

        x1Mask=1:2:length(x);
        x2Mask=2:2:length(x);

        R=getR(x,2,'polar2D');
        
        r=R(:,1);
        theta=R(:,2);
        
    otherwise

end


% get potential parameters
V0=optsPhys.V0;
grav=optsPhys.grav;


VBack = V0*r.^2;

DVBackDr     = 2*V0*r;
DVBackDtheta = zeros(size(theta));

VAdd  = grav*r.*cos(theta).*exp(-0.1*r.^2);

DVAddDr      = grav*cos(theta).*exp(-0.1*r.^2) .* (1 - 2*0.1*r.^2);
DVAddDtheta  = -grav*r.*sin(theta).*exp(-0.1*r.^2);


switch optsPhys.type
    
    case 'stoc'
        
        VAdd  = VAdd(1:2:end);
        VBack = VBack(1:2:end);
        
        DVAdd=zeros(size(x));
        DVBack=zeros(size(x));
        
        DVAdd(x1Mask) = DVAddDr(x1Mask).*cos(theta(x1Mask)) - DVAddDtheta(x1Mask).*sin(theta(x1Mask))./r(x1Mask);
        DVAdd(x2Mask) = DVAddDr(x2Mask).*sin(theta(x2Mask)) + DVAddDtheta(x2Mask).*cos(theta(x2Mask))./r(x2Mask);
        
        r0Mask1 = x1Mask( r(x1Mask)==0 );
        r0Mask2 = x2Mask( r(x2Mask)==0 );
        
        DVAdd(r0Mask1)= DVAddDr(r0Mask1).*cos(theta(r0Mask1));
        DVAdd(r0Mask2)= DVAddDr(r0Mask2).*sin(theta(r0Mask2));

        DVBack(x1Mask) = DVBackDr(x1Mask).*cos(theta(x1Mask)) - DVBackDtheta(x1Mask).*sin(theta(x1Mask))./r(x1Mask);
        DVBack(x2Mask) = DVBackDr(x2Mask).*sin(theta(x2Mask)) + DVBackDtheta(x2Mask).*cos(theta(x2Mask))./r(x2Mask);
        
        DVBack(r==0)=0;
        
    otherwise
        
        DVAdd=[DVAddDr;DVAddDtheta./r];
        DVBack=[DVBackDr;DVBackDtheta./r];
        
end


% add to give total potential and gradient
V1  = VBack + VAdd;
DV1 = DVBack + DVAdd;

end



% figure
% 
% for n=1:length(x1Mask)
%     plot(x1(n),x2(n),'ro')
%     hold on
%     plot([x1(n) x1(n)+DV1(x1Mask(n))],[x2(n) x2(n)+DV1(x2Mask(n))])
% end
% 
% axis 'equal'