function sol = ComputeEquilibrium(this,redo)    

    global dirData
    fprintf('Solving for equilibrium condition...\n');    
    
    %*****************************
    %Initialization
    %*****************************
	rhoLiq_sat    = this.optsPhys.rhoLiq_sat;
	rhoGas_sat    = this.optsPhys.rhoGas_sat;
    kBT           = this.optsPhys.kBT;
    PtsCart     = this.HS.GetCartPts();  
     
    %*****************************
    p         = (this.rho1D_lg-rhoGas_sat)/(rhoLiq_sat-rhoGas_sat);    
    rho_ig    = kron(p,this.rho1D_wl) + kron(1-p,this.rho1D_wg);         
    x_ig      = kBT*log(rho_ig)+this.Vext;
    
    opts             = this.optsNum.PhysArea;
    opts.optsPhys    = this.optsPhys;    
    if(isfield(opts.optsPhys,'gamma'))
        opts.optsPhys    = rmfield(opts.optsPhys,'gamma');
    end
    if(isfield(opts.optsPhys,'Inertial'))
        opts.optsPhys    = rmfield(opts.optsPhys,'Inertial');
    end
    
    if(isfield(opts.optsPhys.V1,'epsilon_w_end'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'epsilon_w_end');
    elseif(isfield(opts.optsPhys.V1,'epsilon_w_Amplitude'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'epsilon_w_Amplitude');
    elseif(isfield(opts.optsPhys.V1,'epsilon_w_max'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'epsilon_w_max');                
    end
    if(isfield(opts.optsPhys.V1,'tau'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'tau');
    end

    opts.maxComp_y2  = this.optsNum.maxComp_y2;
    opts.Comments    = this.configName;
    
    mark             = (PtsCart.y2_kv <= this.optsNum.maxComp_y2);     
            
    %misc = v2struct(mark,Vext,VAdd,Conv,IntMatrFex,x_ig);
    misc.mark = mark;       misc.Vext       = this.Vext; 
    misc.VAdd = this.VAdd;
    misc.x_ig = x_ig;
    misc.Conv = this.Conv;  misc.IntMatrFex = this.IntMatrFex; 
    misc.FexNum = this.optsNum.FexNum;
    
    if(nargin == 1)
        [sol,recEq,paramsEq] = DataStorage('EquilibriumSolutions',...
                            @ComputeEquilibriumCondition,opts,misc); %true      
    else
        [sol,recEq,paramsEq] = DataStorage('EquilibriumSolutions',...
                            @ComputeEquilibriumCondition,opts,misc,redo); %true      
    end

    this.rho_eq      = sol.rho;
    this.FilenameEq  = paramsEq.Filename;
	sol.Filename     = [dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename];
    
    close all;
    
end