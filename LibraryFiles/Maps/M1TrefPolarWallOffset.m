function [y1_kv,y2_kv,J,dH1,dH2] = M1TrefPolarWallOffset(x1,x2,R,angle,ep_angle)

   n = length(x1);
   offset = 0.5;

   c1   = ep_angle*2/(pi*R);
   d    = InvLinearMap(angle,0,pi); 
   %mark = (x1 == -1);
   
   [y1_kv,dy1dx1] =  LinearMap(x1,0,R);  %y1_kv = R*(x1+1)/2;
   
   
    [y2,Diffy2]    = M1Tref(x2,d,c1./(offset+(x1+1).^2));
    y2_kv          = (1+y2)*pi/2;
   
    depdx1         = c1*(-2*(1+x1)./((1+x1).^2+offset).^2); %-c1./(x1.^2+1).^2;
    depddx1        = c1*(8*(1+x1).^2./(offset+(1+x1).^2).^3 - 2./(offset+(1+x1).^2).^2);%2*c1./(x1.^2+1).^3;
      
    dy2dx2         =  Diffy2.dydx*(pi/2);
    dy2ddx2        =  Diffy2.dyddx*(pi/2);
    dy2dx1         =  Diffy2.dyde.*depdx1*(pi/2);
    %dy2dx1(mark)   =  0;
    
    dy2dx1dx2      =  (Diffy2.dydedx.*depdx1)*(pi/2);
   % dy2dx1dx2(mark)=  0;        
    
    dy2ddx1        =  (Diffy2.dydde.*depdx1.^2 + Diffy2.dyde.*depddx1)*(pi/2);
    %x2m            =  x2(mark);         
    %dy2ddx1(mark)  =  ((x2m.^2-1).*(x2m/3-d)./(c1.^2))*(pi/2);
    
    if(nargout >= 3)
        J        = zeros(n,2,2);
        J(:,1,1) = dy1dx1; 
        J(:,1,2) = zeros(n,1);
        J(:,2,1) = dy2dx1;
        J(:,2,2) = dy2dx2;            
    end

    if(nargout >= 4)
        dH1        = zeros(n,2,2);
        dH1(:,1,1) = zeros(n,1); %d2y1dx1dx1; 
        dH1(:,1,2) = zeros(n,1);%d2y1dx1dx2;             
        dH1(:,2,1) = zeros(n,1);%d2y1dx1dx2;     
        dH1(:,2,2) = zeros(n,1);
    end

    if(nargout >= 4)
        dH2        = zeros(n,2,2);
        dH2(:,1,1) = dy2ddx1; %d2y2dx1dx1; 
        dH2(:,1,2) = dy2dx1dx2; %d2y2dx1dx2;             
        dH2(:,2,1) = dy2dx1dx2; %d2y2dx1dx2;     
        dH2(:,2,2) = dy2ddx2; %d2y2dx2dx2;             
    end        
    
end