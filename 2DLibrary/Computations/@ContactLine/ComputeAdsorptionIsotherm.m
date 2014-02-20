 function ComputeAdsorptionIsotherm(this,n)   

    if(isempty(this.HS.Interp))
        InitInterpolation(this);
    end

    [om,rho1Dwg] = Compute1D(this,false,'WG');

    optss                 = this.optsPhys;   
    optss.rho_iguess      = rho1Dwg;
    if((nargin < 2) || isempty(n))
        optss.NContIterations = 200;
    else
        optss.NContIterations = n;
    end

    optss.Dmu = 0;%-0.02;
    
    [rho,ell,mu,OmEx,dmuCheck] = FMT_1DContinuation(this.HS,...
        this.IntMatrFex,optss,this.optsNum.FexNum,this.Conv,'mu','load');%,'Movie');            

    this.AdsorptionIsotherm_FT       = ell; 
    this.AdsorptionIsotherm_Mu       = mu;
    this.AdsorptionIsotherm_rho      = rho;
    this.AdsorptionIsotherm_OmEx     = OmEx;
    this.AdsorptionIsotherm_dmuCheck = dmuCheck;                       

end 