function [rho_ic1D,postParms] = FMT_1D_Dyn(HS,IntMatrFex_2D,optsPhys,FexNum,Conv)
%external background potential =0
    
    close all;      

    %************************************************
    %***************  Initialization ****************
    %************************************************
    kBT       = optsPhys.kBT;           
    R         = optsPhys.sigmaS/2;
    fBulk     = str2func(['FexBulk_',FexNum.Fex]);  
    
    if(isfield(optsPhys,'Dmu') && isfield(optsPhys,'mu_sat'))
        mu        = optsPhys.mu_sat + optsPhys.Dmu;
    elseif(isfield(optsPhys,'eta'))
        eta       = optsPhys.eta;
        rhoBulk   = eta*6/pi;
        mu        = kBT*log(rhoBulk) + fBulk(rhoBulk,kBT);
    end
    
    getFex = str2func(['Fex_',FexNum.Fex]);
    
    markComp  = (HS.Pts.y1_kv==inf);        
    nSpecies  = 1;
    Pts       = HS.Pts;
    N         = sum(markComp);
    %************************************************
    %****************  Preprocess  ******************
    %************************************************
        
    close all;
    IntMatrFex_1D  = Get1DMatrices(IntMatrFex_2D,HS);    
    IntMatrFex     = IntMatrFex_1D;
    
    if(isempty(Conv))
        Conv = zeros(size(IntMatrFex));
    else
        Conv = Conv(markComp,markComp);
    end    
           
    y1S     = Pts.y1_kv(markComp); 
    y2S     = Pts.y2_kv(markComp);
    
    ptsCart = HS.GetCartPts(y1S,y2S);
    VAdd    = getVAdd(ptsCart.y1_kv,ptsCart.y2_kv,0,optsPhys.V1);    
    
	[h1s,h2s,Int_1D]     = HS.ComputeIntegrationVector();
    [h1s,h2s,Int_1D_AD]  = HS.AD.ComputeIntegrationVector();
    
    Diff.Dy  = HS.Diff.Dy2(markComp,markComp);
    Diff.DDy = HS.Diff.DDy2(markComp,markComp);    

    %****************************************************************
    %**************** Solve for equilibrium 1D condition   **********
    %****************************************************************        
    tic
    if(isfield(optsPhys,'rho_iguess'))
        y0        = kBT*log(optsPhys.rho_iguess)*ones(length(Pts.y1_kv),nSpecies);
    else
        y0        = zeros(length(Pts.y1_kv),nSpecies);%getInitialGuess(VAdd0);                
    end       
    y0 = y0(markComp);
    
    %PlotRosenfeldFMT_AverageDensities(HS,IntMatrFex(1),ones(size(y0)));                       
    fsolveOpts=optimset('MaxFunEvals',2000000,'MaxIter',200000,'Display','off');    
    [x_ic_1D,h1,flag] = fsolve(@f,y0,fsolveOpts);     
    if(flag ~= 1)
        cprintf('red','Error in fsolve, FMT_1D_Interface');
    end
    rho_ic1D  = exp(x_ic_1D/kBT); 
    
    %**************************************
    %**************** Dynamics   **********
    %**************************************    
     if(optsPhys.Inertial)
         mMX     = ones(N,1);    
         mMV     = ones(N,1);
         mMV(1)   = 0;
         mMX(end) = 0;
         mM       = [mMX;mMV];                   
         opts = odeset('RelTol',10^-10,'AbsTol',10^-10,'Mass',diag(mM));
         [outTimes,X_t] =  ode15s(@dx_dtInertia,optsPhys.plotTimes,...
                             [x_ic_1D;zeros(N,1)],opts);
                         
        for i=1:length(outTimes)
            rho_ic1D_t  = exp(X_t(i,1:end/2)'/kBT);
            v_t         = X_t(i,1+end/2:end)';
            
            subplot(1,2,1);
            HS.do1DPlotNormal(rho_ic1D_t);
            title(['t = ',num2str(outTimes(i))]);
            
            subplot(1,2,2);
            hold off;
            HS.do1DPlotNormal(v_t);
            pause(0.4);
        end                         
     else        
        mM      = ones(N,1);    
        mM(1)   = 0;  
        mM(end) = 0;        
        opts = odeset('RelTol',10^-10,'AbsTol',10^-10,'Mass',diag(mM));
        [outTimes,X_t] =  ode15s(@dx_dt,optsPhys.plotTimes,x_ic_1D,opts);
        
        for i=1:length(outTimes)
            rho_ic1D_t  = exp(X_t(i,:)'/kBT);
            HS.do1DPlotNormal(rho_ic1D_t);
            title(['t = ',num2str(outTimes(i))]);
            pause(0.4);
        end
     end        
    
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************             
    function y = f(x)
        %solves for T*log*rho                        
        y            = GetExcessChemPotential(x,0,mu);         
        y            = y(:);
    end
    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp(x/kBT);                
        mu_s  = getFex(rho_s,IntMatrFex,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end
        VAdd    = getVAdd(ptsCart.y1_kv,ptsCart.y2_kv,t,optsPhys.V1);            
        mu_s = mu_s + x + Conv*rho_s + VAdd;
    end

    function dxdt = dx_dt(t,x)              
        
        mu_s     = GetExcessChemPotential(x,t,mu);        
        mu_s(end) = 0;
        h_s      = Diff.Dy*x;
                
        dxdt     = kBT*Diff.DDy*mu_s + h_s.*(Diff.Dy*mu_s);
                
        %Boundary Conditions at infinity
        flux_dir    = Diff.Dy*mu_s;
        dxdt(1)     = flux_dir(1);
        dxdt(end)   = x(end)-x_ic_1D(end);

        dxdt = dxdt(:);
    end 
    function dydt = dx_dtInertia(t,y)
        
        x = y(1:end/2);
        v = y(1+end/2:end);        
        
        mu_s     = GetExcessChemPotential(x,t,mu);        
        h_s      = Diff.Dy*x;
                
        dxdt     = - kBT*(Diff.Dy*v) - v.*h_s;
        dvdt     = - v.*(Diff.Dy*v) - optsPhys.gamma*v - Diff.Dy*mu_s;                
                
        %Boundary Conditions at infinity        
        dvdt(1)     = v(1);
        dxdt(end)   = x(end)-x_ic_1D(end);

        dydt = [dxdt;dvdt];
    end 
end    