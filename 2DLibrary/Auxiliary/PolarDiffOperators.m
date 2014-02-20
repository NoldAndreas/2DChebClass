function DiffPolar = PolarDiffOperators(Diff,Pts)        
    %Please note: entries including 1/r are set to zero at infinity! 
    % This is for convenience of using the term d/dr in grad up to the
    % origin. 

    N1 = length(Pts.x1);  N2 = length(Pts.x2);
    
    %First order operators
    DR      = Diff.Dy1;       
	Dtheta  = Diff.Dy2;     
    
    Diag_r1  = diag(1./Pts.y1_kv);
    Diag_r1(Diag_r1 == inf) = 0;
    
    grad           = [DR ; Diag_r1*Dtheta];  %gradient of a scalar
    DRDtheta       = Diff.Dy1Dy2;
    %at origin, 1/r*(d/dtheta) = d^2/(dtheta dr), given that d/dtheta = 0
    grad([false(N1*N2,1);Pts.y1_kv==0],:) = DRDtheta(Pts.y1_kv==0,:);
    
    div            = [DR  + Diag_r1 , Diag_r1*Dtheta]; %divergence of a vector
    div(Pts.y1_kv == 0,[true(N1*N2,1);false(N1*N2,1)]) = 2*DR(Pts.y1_kv==0,:);
    div(Pts.y1_kv == 0,[false(N1*N2,1);true(N1*N2,1)]) = DRDtheta(Pts.y1_kv==0,:);    
    
    
    if(~isfield(Diff,'DDy1'))
        DiffPolar = struct('Dy1',Diff.Dy1,'Dy2',Diff.Dy2,'grad',grad,'div',div);                
        return;
    end
    
    %Second Order
    DDR     = Diff.DDy1;
    DDtheta = Diff.DDy2;
    
    Diag_r2  = Diag_r1.^2;
        
    % in the following, 'Vector' will denote (u,v)', where u is a (e.g.)
    % velocity in radial direction, v in angular direction
    rdiag   = diag(Pts.y1_kv);
    %rdiag(rdiag == inf) = 0;
    %LapR    = rdiag*DDR + DR + Diag_r1*DDtheta;  %R * Laplace of a Scalar                       
    Lap     = DDR + Diag_r1*DR + Diag_r2*DDtheta;  %Laplace of a Scalar                       
    %at origin, Lap = 2d^2/dr^2 + 1/r, given that for r= 0,
    % (1) d^2/dtheta^2 = 0
    % (2) d/dr + d^3/
    help_ddrddth = DDtheta*DDR;
    Lap(Pts.y1_kv == 0,:) = 2*DDR(Pts.y1_kv==0,:) + 1/2*help_ddrddth(Pts.y1_kv==0,:);
    
    
    LapVec  = [Lap - Diag_r2     , -2*Diag_r2* Dtheta; %Laplace of a Vector
                   2*Diag_r2*Dtheta  , Lap - Diag_r2];    
               
               
    if(~isfield(Diff,'DDDy1'))
        DiffPolar = struct('Dy1',Diff.Dy1,'DDy1',Diff.DDy1,...
                      'Dy2',Diff.Dy2,'DDy2',Diff.DDy2,...
                      'Dy1Dy2',Diff.Dy1Dy2,...
                      'Lap',Lap,'LapVec',LapVec,'grad',grad,'div',div);        
        return;
    end
    
    %Third order
    DDDR    = Diff.DDDy1;
    DDDtheta= Diff.DDDy2;                    
        
    DRDDtheta = Diff.Dy1DDy2;
    DDRDtheta = Diff.DDy1Dy2;
    
        
    Diag_r3   = Diag_r1.^3;       
        
   
    % The following operators are not defined successively for reasons of
    %accuracy:
    
    %gradient of the Laplace of a Scalar     
    gradLap     = [DDDR - Diag_r2*DR + Diag_r1*DDR - 2*Diag_r3*DDtheta + Diag_r2*DRDDtheta ;...
                        Diag_r1*DDRDtheta + Diag_r2*DRDtheta + Diag_r3*DDDtheta]; 
   
    %gradient of the Divergence of a Vector
    gradDivVec  = [DDR - Diag_r2 + Diag_r1*DR , - Diag_r2*Dtheta + Diag_r1*DRDtheta;...
                           Diag_r1*DRDtheta+Diag_r2*Dtheta , Diag_r2*DDtheta];
                   
    g1          = DDDR + Diag_r3 + Diag_r3*DDtheta - Diag_r2*DR + 2*Diag_r1*DDR + Diag_r2*DRDDtheta;
    g2          = Diag_r3*(Dtheta + DDDtheta) + Diag_r1*DDRDtheta - Diag_r2*DRDtheta;
    LapDiv      = [g1 g2];
    
    DiffPolar   = struct('Dy1',Diff.Dy1,'DDy1',Diff.DDy1,'DDDy1',Diff.DDDy1,...
                      'Dy2',Diff.Dy2,'DDy2',Diff.DDy2,'DDDy2',Diff.DDDy2,...
                      'Dy1Dy2',Diff.Dy1Dy2,'DDy1Dy2',Diff.DDy1Dy2,'Dy1DDy2',Diff.Dy1DDy2,...
                      'Lap',Lap,'LapVec',LapVec,'grad',grad,'div',div,'gradLap',gradLap,...
                      'gradDivVec',gradDivVec,'LapDiv',LapDiv);
end