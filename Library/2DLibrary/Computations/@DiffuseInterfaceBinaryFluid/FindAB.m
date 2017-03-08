function vec_a = FindAB(this,phi,mu,deltaX,a_ig)
    
    Ind     = this.IDC.Ind;
    M       = this.IDC.M;        
    F       = false(M,1);    
    IBB     = repmat(Ind.bound,2,1);
    PtsCart = this.IDC.GetCartPts();    
    
    if(isempty(this.theta))
        theta = pi/2;
    else
        theta = this.theta;
    end            
    
    if(nargin == 1)
        if(isempty(this.phi))
            phi     = InitialGuessRho(this,theta);                
            mu       = zeros(M,1);
            p       = zeros(M,1);
        else
            phi = this.phi;            
            mu   = this.mu;
            p   = this.p;
        end

        if(isempty(this.deltaX))       
            deltaX  = 0;        
        else
            deltaX  = this.deltaX;
        end

        if(~isempty(this.a))
            a_ig    = this.a;
        else
            a_ig    = [0;0];
        end            
    end
    
  %  CheckHalfLineIntegral();
    initialError = f(a_ig);
    disp(['Initial value: ',num2str(a_ig(1)),',',num2str(a_ig(2)),'.']);
    disp(['Initial error: ',num2str(initialError(1)),',',num2str(initialError(2)),'.']);
    vec_a = fsolve(@f,a_ig)
        
    uv(IBB)     = GetVelBC(this,zeros(2*M,1),vec_a,deltaX,theta);
    uv_org(IBB) = GetVelBC(this,zeros(2*M,1),[0;0],deltaX,theta);
    
    plot(PtsCart.y1_kv(Ind.top),-uv([Ind.top;F]),'b'); hold on;
    plot(PtsCart.y1_kv(Ind.top),-uv([F;Ind.top]),'r');
    
    plot(PtsCart.y1_kv(Ind.top),-uv_org([Ind.top;F]),'b--'); hold on;
    plot(PtsCart.y1_kv(Ind.top),-uv_org([F;Ind.top]),'r--');
    xlim([-40,40]);
    
    function v = f(x)
        a       = x(1:2);        
        uv_     = zeros(2*M,1);        
        
        uv_(IBB)   = GetVelBC(this,uv_,a,deltaX,theta);                
        v_SeppAdd  = GetSeppecherConditions(this,-uv_,phi,mu,zeros(M,1),deltaX,theta);

        v = v_SeppAdd(1:2);        
    end

    function CheckHalfLineIntegral()
        
        N     = 50;
        y2Max = 10;
        y10   = y2Max/tan(theta);
        
        optsPlot = struct('yMin',y10,'yMax',y10+20,'N',200);
        
        opts = struct('L',10,'N',N,'yMin',y10);
        HIS = HalfInfSpectralLine(opts);
        HIS.ComputeAll(optsPlot);
        Int = HIS.Int;        
        
        Pts.y1_kv = HIS.Pts.y;
        Pts.y2_kv = y2Max*ones(size(Pts.y1_kv));
        
        u_flow    = GetSeppecherSolutionCart([Pts.y1_kv,...
                                              Pts.y2_kv],1,0,0,theta);                  
    
        disp(['Error of right part: ',num2str(Int*u_flow(1+end/2:end)+y2Max)]);
        HIS.plot(u_flow(1+end/2:end)); hold on;
        
        u_flow    = GetSeppecherSolutionCart([2*y10-Pts.y1_kv,...
                                              Pts.y2_kv],1,0,0,theta);                  
    
        disp(['Error of left part: ',num2str(Int*u_flow(1+end/2:end)-y2Max)]);
       
        opts     = struct('L',7,'N',N);
        optsPlot = struct('yMin',-100,'yMax',100,'N',200);
        ISL      = InfSpectralLineQuad(opts);
        
        ISL.ComputeAll(optsPlot);
        Int  = ISL.Int;
        
        Pts.y1_kv = ISL.Pts.y;
        Pts.y2_kv = y2Max*ones(size(Pts.y1_kv));
        
        u_flow    = GetSeppecherSolutionCart([Pts.y1_kv,...
                                              Pts.y2_kv],1,0,0,theta);                  
    
        disp(['Error of full part: ',num2str(Int*u_flow(1+end/2:end))]);
        v = u_flow(1+end/2:end);
        ISL.plot(v); hold on;
        plot([y10,y10],[min(v) max(v)]);
    end
end