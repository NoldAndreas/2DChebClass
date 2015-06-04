 function ComputeEquilibrium(this,rho_ig,optsIn,miscIn)
    
    fprintf('Solving for equilibrium condition...\n');    
    
    %*****************************
    %Initialization
    %*****************************	    
    if((nargin == 1) || isempty(rho_ig))
        x_ig         = getInitialGuess(this);
    else
        x_ig         = getInitialGuess(this,rho_ig);
    end
    if(nargin>2)
        opts = optsIn;
    end
    if(nargin>3)
        misc = miscIn;
    end
    
    opts.PhysArea    = this.optsNum.PhysArea;
    opts.optsPhys    = this.optsPhys;    
    
    % Clean optsPhys from unnecessary information
    if(isfield(opts.optsPhys,'Inertial'))
        opts.optsPhys    = rmfield(opts.optsPhys,'Inertial');
    end    
    if(isfield(opts.optsPhys,'HI'))
        opts.optsPhys    = rmfield(opts.optsPhys,'HI');
    end
    if(isfield(opts.optsPhys,'nParticles'))
        opts.optsPhys    = rmfield(opts.optsPhys,'nParticles');
    end
    
    %opts.maxComp_y2  = this.optsNum.maxComp_y2;
%    opts.Comments    = this.configName;    
    %mark             = (PtsCart.y2_kv <= this.optsNum.maxComp_y2);     
            
    %misc = v2struct(mark,Vext,VAdd,Conv,IntMatrFex,x_ig);
    %misc.mark = mark;           
    misc.Vext       = this.Vext; 
    misc.VAdd       = this.VAdd;
    misc.x_ig       = x_ig;
    misc.Conv       = this.IntMatrV2;  
    misc.IntMatrFex = this.IntMatrFex;
    misc.FexNum     = this.optsNum.FexNum;
    misc.Int        = this.IDC.Int;
    
    opts.Comments = this.configName;
    
    ignoreList = {'optsPhys_p','optsPhys_mu_sat','optsPhys_rhoLiq_sat',...
                  'optsPhys_rhoGas_sat',...
                  'optsPhys_gamma','optsPhys_gammaS','optsPhys_D0'...
                  'optsPhys_tMax','D0S',...
                  'optsPhys_V1_epsilon_w_end',...
                  'optsPhys_V1_epsilon_w_Amplitude',...
                  'optsPhys_V1_epsilon_w_max',...
                  'optsPhys_V1_tau',...
                  'optsPhys_viscosity'};
	[sol,recEq,paramsEq] = DataStorage('Equilibrium',...
                            @ComputeEquilibriumCondition,opts,misc,[],ignoreList);
         
    this.x_eq   = sol.x;
    this.mu     = sol.mu;
    if(isfield(paramsEq,'Filename'))
        this.FilenameEq  = paramsEq.Filename;
    end
    
end