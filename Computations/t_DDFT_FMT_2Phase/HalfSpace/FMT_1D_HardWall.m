function rho_ic1D = FMT_1D_HardWall(HS,IntMatrFex_2D,optsPhys,optsNum)
%************************************************************************* 
% define
%   mu_s_i := kBT*log(rho_i) + sum( int(rho_j(r')*Phi2D(r-r'),dr') , 
%               j = 1..N) +mu_HS(rho_i) + V_ext_i - mu_i
% Equilibrium, (solve for each species i)
%  (EQ i,1) mu_s_i     = 0
%  (EQ i,2) int(rho_i) = NoParticles_i
% Dynamics: 
%*************************************************************************       
    %********************************************
    %**************** Initialization   **********
    %********************************************
    
%    saveFigs  = true;
%    eta       = optsPhys.eta;    
    
    markComp  = (HS.Pts.y1_kv==inf);    
	Pts                       = HS.Pts;
    PtsCart                   = HS.GetCartPts();
    Interp1D                  = HS.ComputeInterpolationMatrix(1,(-1:0.01:0.7)',true,true);        
    Interp1D.InterPol         = Interp1D.InterPol(:,markComp);
    Interp1D.ptsCart          = HS.GetCartPts(Interp1D.pts1,Interp1D.pts2);
         
    subPts.y2_kv              = (0.:0.01:3.5)';
    subPts.y1_kv              = inf*ones(size(subPts.y2_kv));           
    IP                        = HS.AD.SubShapePts(subPts);
    Interp1D_AD.InterPol      = IP(:,HS.AD.Pts.y1_kv == inf);
    Interp1D_AD.pts1          = subPts.y1_kv; 
    Interp1D_AD.pts2          = subPts.y2_kv;
    
%	subPtsAAD.y2_kv           = (0.5:0.01:3.5)';
%    subPtsAAD.y1_kv           = inf*ones(size(subPtsAAD.y2_kv));           
%    IP                        = HS.SubShapePts(subPtsAAD);
%    Interp1D_AAD.InterPol     = IP(:,markComp);
%    Interp1D_AAD.pts1         = subPtsAAD.y1_kv; 
%    Interp1D_AAD.pts2         = subPtsAAD.y2_kv;

	%***********************************
    %**************** Compute **********
    %***********************************
    rho_ic1D = FMT_1D(HS,IntMatrFex_2D,optsPhys,optsNum.FexNum,[],true);
    

end