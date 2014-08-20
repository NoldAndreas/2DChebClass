function [f,y1] = PostProcess(this,y1Int,y2Int)
        
    InitAnalysisGrid(this,y1Int,y2Int);
    ComputeInterfaceContour(this);
    
    ComputeAdsorptionIsotherm(this);	    
    disjoingPressure1DCheck();    
    
    PostProcess_2DDisjoiningPressure(this);   
    SumRule_DisjoiningPotential(this);
    
    %f = f2;
    %f = fB;            
%    this.grandPot      = f;   
%    this.grandPot2     = f2;           
    
  %  this.wIntLoc       = wIntLoc;
   % this.wIntHS        = wIntHS;
    %Compute Film thickness as contour line       
    
    function disjoingPressure1DCheck()
        rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
        rhoGas_sat     = this.optsPhys.rhoGas_sat;
        kBT            = this.optsPhys.kBT;
        mu_sat         = this.optsPhys.mu_sat;    
        Dmu            = this.optsPhys.Dmu;    
        R              = this.optsPhys.sigmaS/2; 
        getFex         = str2func(['Fex_',this.optsNum.FexNum.Fex]);    

        optsPhys = this.optsPhys;        

        rho_Bulk       = rhoGas_sat;
        mu             = mu_sat+Dmu;

        rho            = this.rho_eq;    

        %******************************************

        %Get Bulk value and compare:        
        [h1s,f_hs_Bulk,h2s,h3s] = FexBulk_FMTRosenfeld_3DFluid(rho_Bulk,kBT);
        
        Phi_r             = str2func(optsPhys.V2.V2DV2);        
        [h1s,h2s,alpha]   = Phi_r(0);
        f_attr_Bulk       = alpha*rho_Bulk^2;
        f_id_Bulk         = kBT*rho_Bulk.*(log(rho_Bulk)-1);
        f_Vmu_Bulk        = -mu*rho_Bulk;

        floc_Bulk         = f_id_Bulk + f_attr_Bulk + f_Vmu_Bulk;

        %Compute excess grand potential        
        f_id           = kBT*rho.*(log(rho)-1);
        [h1s,h2s,f_hs] = getFex(rho,this.IntMatrFex,kBT,R);
        f_attr         = 0.5*rho.*(this.Conv*rho);
                
        f_Vmu          = rho.*(this.VAdd- mu);        

        f_loc  = f_id + f_attr + f_Vmu;        

        this.f_loc = f_loc;
        this.f_hs  = f_hs;

        %Initialize Integration Boxes which allow for integration normal to
        %wall
        y2Max = this.optsNum.maxComp_y2 + 10;

        %y1 = this.y1;
    %    f             = zeros(size(y1));
    %    f2            = zeros(size(y1));

      %  wIntLoc       = zeros(length(y1),length(rho));
     %   wIntHS        = zeros(length(y1),length(f_hs));

        PtsCart       = this.HS.GetCartPts();   
        [h1,dVAdd]    = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,0,this.optsPhys.V1);


        [I,w,weights,IP,pts] = this.HS.doIntFLine([inf inf],[0.5 y2Max],[],'CHEB');
        [h1,dVadd_i]         = getVAdd(pts.y1_kv,pts.y2_kv,0,this.optsPhys.V1);
        IP                   = IP(:,this.HS.Pts.y1_kv == inf);
        fB_AI                = zeros(size(this.AdsorptionIsotherm_FT));
        for i = 1:length(this.AdsorptionIsotherm_FT)
            optss              = this.optsPhys;
            optss.Dmu          = this.AdsorptionIsotherm_Mu(i) - mu_sat;
            optss.rho_iguess   = 'WL';
            [rho_wl,postParms] = FMT_1D(this.HS,this.IntMatrFex,optss,this.optsNum.FexNum,this.Conv,false);

            rho_i    = this.AdsorptionIsotherm_rho(i,:);
            fB_AI(i) = - weights*(dVadd_i.dy2.*(IP*(rho_i'-rho_wl))) + ...
                             + this.optsPhys.kBT*(rho_i(1) - rho_wl(1));                
            %dmuCheck(i) = -(kBT*(rho_i(1) - rho_wl(1)) - Int_1D*( (rho_i - rho_wl).*dVAdd.dy2))/(rhoLiq_eq-rhoGas_eq);
        end

        this.disjoiningPressureCheck = fB_AI;    
    end
end