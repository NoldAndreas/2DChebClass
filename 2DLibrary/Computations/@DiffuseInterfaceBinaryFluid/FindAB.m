function FindAB(this)
    
    Ind     = this.IDC.Ind;
    M       = this.IDC.M;    
    T       = true(M,1);
    F       = false(M,1);
    ITT     = repmat(Ind.top,2,1);
    IBB     = repmat(Ind.bound,2,1);
    
    deltaX  = 0;
    theta   = 1.4;
    phi     = InitialGuessRho(this);
    G       = zeros(M,1);
    p       = zeros(M,1);
    PtsCart = this.IDC.GetCartPts;
        
    vec = fsolve(@f,[0;0]) 
    
    uv(IBB)     = GetVelBC(this,uv,a,deltaX,theta);
    uv_org(IBB) = GetVelBC(this,uv,[0;0],deltaX,theta);
    
    plot(PtsCart.y1_kv(Ind.top),-uv([Ind.top;F]),'b'); hold on;
    plot(PtsCart.y1_kv(Ind.top),-uv([F;Ind.top]),'r');
    
    plot(PtsCart.y1_kv(Ind.top),-uv_org([Ind.top;F]),'b--'); hold on;
    plot(PtsCart.y1_kv(Ind.top),-uv_org([F;Ind.top]),'r--');
    xlim([-10,10]);
    function v = f(x)

        a       = x(1:2);        
        uv      = zeros(2*M,1);        
        
        uv(IBB)    = GetVelBC(this,uv,a,deltaX,theta);                
        v_SeppAdd  = GetSeppecherConditions(this,-uv,phi,G,p,deltaX,theta);

        v = v_SeppAdd(1:2);        
    end
end