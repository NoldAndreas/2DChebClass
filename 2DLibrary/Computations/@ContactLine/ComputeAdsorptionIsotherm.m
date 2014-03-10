 function ComputeAdsorptionIsotherm(this,n)   

    if(isempty(this.HS.Interp))
        InitInterpolation(this);
    end

    [om,rho1Dwg] = Compute1D(this,false,'WG');

        optss              = this.optsPhys;   
        optss.rho_iguess   = rho1Dwg;
        optss.Dmu          = 0;
        
        if((nargin < 2) || isempty(n))
            optss.NContIterations = 200;
        elseif(isnumeric(n))
            optss.NContIterations = n;        
        end

        if(strcmp(n,'load'))
            [rho,ell,mu,OmEx,dmuCheck,pts] = FMT_1DContinuation(this.HS,...
                this.IntMatrFex,optss,this.optsNum.FexNum,this.Conv,'mu','load');%,'Movie');            
        else
            [rho,ell,mu,OmEx,dmuCheck,pts] = FMT_1DContinuation(this.HS,...
                this.IntMatrFex,optss,this.optsNum.FexNum,this.Conv,'mu');%,'load');%,'Movie');            
        end
    
    this.AdsorptionIsotherm_FT       = ell; 
    this.AdsorptionIsotherm_Mu       = mu;
    this.AdsorptionIsotherm_rho      = rho;
    this.AdsorptionIsotherm_OmEx     = OmEx;
    this.AdsorptionIsotherm_dmuCheck = dmuCheck;  
    this.AdsorptionIsotherm_Pts      = pts;

end 