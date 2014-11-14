function FindAB(this)
    
    Ind     = this.IDC.Ind;
    M       = this.IDC.M;    
    T       = true(M,1);
    F       = false(M,1);
    ITT     = repmat(Ind.top,2,1);
    IBB     = repmat(Ind.bound,2,1);
    
    deltaX  = 0;
    theta   = 100*pi/180;
    phi     = InitialGuessRho(this);
    G       = zeros(M,1);
    p       = zeros(M,1);
    PtsCart = this.IDC.GetCartPts;
        
  %  CheckHalfLineIntegral();
    
    vec = fsolve(@f,[0;0])    
        
    uv(IBB)     = GetVelBC(this,zeros(2*M,1),a,deltaX,theta);
    uv_org(IBB) = GetVelBC(this,zeros(2*M,1),[0;0],deltaX,theta);
    
    plot(PtsCart.y1_kv(Ind.top),-uv([Ind.top;F]),'b'); hold on;
    plot(PtsCart.y1_kv(Ind.top),-uv([F;Ind.top]),'r');
    
    plot(PtsCart.y1_kv(Ind.top),-uv_org([Ind.top;F]),'b--'); hold on;
    plot(PtsCart.y1_kv(Ind.top),-uv_org([F;Ind.top]),'r--');
    xlim([-40,40]);
    function v = f(x)

        a        = x(1:2);        
        uv_      = zeros(2*M,1);        
        
        uv_(IBB)   = GetVelBC(this,uv_,a,deltaX,theta);                
        v_SeppAdd  = GetSeppecherConditions(this,-uv_,phi,G,p,deltaX,theta);

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
        HIS.doPlots(u_flow(1+end/2:end)); hold on;
        
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
        ISL.doPlots(v); hold on;
        plot([y10,y10],[min(v) max(v)]);
    end
end