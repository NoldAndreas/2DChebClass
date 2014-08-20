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
    if(isfield(opts.optsPhys,'gamma'))
        opts.optsPhys    = rmfield(opts.optsPhys,'gamma');
    end
    if(isfield(opts.optsPhys,'gammaS'))
        opts.optsPhys    = rmfield(opts.optsPhys,'gammaS');
    end
    if(isfield(opts.optsPhys,'D0'))
        opts.optsPhys    = rmfield(opts.optsPhys,'D0');
    end
    if(isfield(opts.optsPhys,'D0S'))
        opts.optsPhys    = rmfield(opts.optsPhys,'D0S');
    end
    if(isfield(opts.optsPhys,'Inertial'))
        opts.optsPhys    = rmfield(opts.optsPhys,'Inertial');
    end    
    if(isfield(opts.optsPhys,'HI'))
        opts.optsPhys    = rmfield(opts.optsPhys,'HI');
    end
    if(isfield(opts.optsPhys,'tMax'))
        opts.optsPhys    = rmfield(opts.optsPhys,'tMax');
    end
    if(isfield(opts.optsPhys,'nParticles'))
        opts.optsPhys    = rmfield(opts.optsPhys,'nParticles');
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
    
	[sol,recEq,paramsEq] = DataStorage('EquilibriumSolutions',...
                            @ComputeEquilibriumCondition,opts,misc);
         
    this.x_eq   = sol.x;
    this.mu     = sol.mu;
    this.FilenameEq  = paramsEq.Filename;
    
end