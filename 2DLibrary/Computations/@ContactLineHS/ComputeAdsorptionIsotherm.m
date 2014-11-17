 function ComputeAdsorptionIsotherm(this,n,drying)           
	optss              = this.optsPhys;   
    
    if((nargin > 2))
            optss.drying = drying;        
	else
        if(IsDrying(this))
            optss.drying = 'drying';        
        else
            optss.drying = 'wetting';
        end 
    end
    
    optss.Dmu          = 0;

    if((nargin < 2) || isempty(n))
        optss.NContIterations = 200;
    elseif(isnumeric(n))
        optss.NContIterations = n;        
    end

%    if(strcmp(n,'load'))                            
%       [rho,ell,mu,OmEx,dmuCheck,pts] = FMT_1DContinuation(this.IDC,...
%            this.IntMatrFex,optss,this.optsNum.FexNum,this.IntMatrV2.Conv,'mu','load');%,'Movie');            
%    else
        [rho,ell,mu,OmEx,dmuCheck,pts] = FMT_1DContinuation(this.IDC,...
            this.IntMatrFex,optss,this.optsNum.FexNum,this.IntMatrV2.Conv,'mu');%,'load');%,'Movie');            

%    end

    this.AdsorptionIsotherm = struct('FT',ell,...
                                 'mu',mu,...
                                 'rho',rho,...
                                 'OmEx',OmEx,...
                                 'dmuCheck',dmuCheck,...
                                 'Pts',pts);
    
end 