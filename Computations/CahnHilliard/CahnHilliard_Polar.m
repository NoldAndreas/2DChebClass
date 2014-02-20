function data = CahnHilliard_Polar(optsPhys,optsNum)
%************************************************************************* 
%data = CahnHilliard_Polar(optsPhys,optsNum)
%
%            (Eq)  0 = 1/Cn * d(rho*f)/drho -Cn*Lap(rho) - mu + V_ext
%                    = (rho*(1-rho)*(1-2rho))/Cn -Cn*Lap(rho) - mu + V_ext
% with   (BC Wall) 0 = Cn*normal*grad(rho_s) + cos(theta)*(rho-1)*rho
%(BC outer Radius) 0 = normal*grad(rho_s);
%
%*************************************************************************   
    if(nargin == 0)    
        
        Phys_Area = struct('y1Min',0,'y1Max',10,'N',[30;20],...
                       'y2Min',-pi/2,'y2Max',pi/2);
        
        Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,...
                       'y2Min',-pi/2,'y2Max',pi/2,'N2',100);

        Sub_Area = struct('y1Min',-5,'y1Max',2,'N',[50,50],...
                           'y2Min',-2,'y2Max',2);        

        optsNum = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'SubArea',Sub_Area,...
                         'DDFTCode','CahnHilliard_Polar',...
                         'plotTimes',0:0.1:5);                            

        optsPhys = struct('V1DV1','Vext_Cart_5',...
                          'V0',0.0,'epsilon_w',0,'epsilon_w_end',0,'y10',2,'y20',0,'tau',1,...
                          'theta',-pi/3,'Cn',0.5,...
                         'kBT',0.7,...
                         'nParticlesS',pi*100/8);
    end

    disp(['** ',optsNum.DDFTCode,' **']);
 
    %************************************************
    %***************  Initialization ****************
    %************************************************
    close all;
    
    kBT         = optsPhys.kBT; 
    nParticles  = optsPhys.nParticlesS;            
    plotTimes   = optsNum.plotTimes;
    Cn          = optsPhys.Cn;    
    theta       = optsPhys.theta;
        
    PhysArea    = optsNum.PhysArea;
    N1 = PhysArea.N(1); N2 = PhysArea.N(2);                                       
                  
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    tic
    WDG              = Wedge(PhysArea);
    [Pts,Diff,Int,Ind,Interp] = WDG.ComputeAll(optsNum.PlotArea);        
    
    %subBox            = Box(optsNum.SubArea);    
    %IP                = adisc.SubShapePts(subBox.Pts);
    %Int_SubOnFull     = subBox.ComputeIntegrationVector()*IP;

%    [Path,InterpPath,Int_of_path,Int_SubOnFull] = SubSpace(SubArea,...
%                @SpectralSpectral_Interpolation,Pts,Maps,'normal');                   
    
%    t_preprocess = toc;
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    tic
    [Vext,Vext_grad]  = getVBackDVBack(Pts.y1_kv,Pts.y2_kv,optsPhys);
    mu                = 0;
    rho_ic            = fsolve(@f_eq,InitialGuess());
    
    subplot(2,2,1); WDG.doPlots(rho_ic,'SC');
    [fE,dfE]          = f(rho_ic);
    subplot(2,2,2); WDG.doPlots(-(rho_ic.^2).*dfE/Cn,'SC');
    subplot(2,2,3); WDG.doPlots(rho_ic.*fE/Cn-mu*rho_ic,'SC');    
%    t_eqSol = toc;
%     %****************************************************************
%     %****************  Compute time-dependent solution   ************
%     %****************************************************************
%     tic    
%     mM             = ones(N1*N2,1);
%     mM(Ind.bound)  = 0;
%     opts           = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag([mM]));
%     %[outTimes,X_t] =  ode15s(@dx_dt,plotTimes,[0;rho_ic],opts);        
%     [outTimes,X_t] =  ode15s(@dx_dt,plotTimes,rho_ic,opts);        
%     t_dynSol = toc;
%     %************************************************
%     %****************  Postprocess  ****************
%     %************************************************        
% 
%     accFlux   = X_t(:,1);
%     X_t       = X_t(:,2:end)';
%     rho_t     = exp((X_t-Vext*ones(1,length(outTimes)))/kBT);
%     
%     flux_t    = zeros(2*N1*N2,length(plotTimes));
%     for i = 1:length(plotTimes)
%         flux_t(:,i) = GetFlux(X_t(:,i),plotTimes(i));
%     end
%     
%     Subspace = struct('subArea',subdisc);%Path,InterpPath,Int_of_path,Int_SubOnFull,accFlux);
%     data     = v2struct(Conv,X_t,rho_t,Subspace,...
%                         mu,flux_t,t_preprocess,t_eqSol,t_dynSol);
%     data.shape = adisc;    
%     
%     if(~isfield(optsNum,'savefileDDFT'))
%         SaveToFile(optsNum.DDFTCode,v2struct(data,optsPhys,optsNum),getResultsPath());
%     end                    
%                     
% 	display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess)]);
%     display(['Equilibrium Sol., Computation time (sec): ', num2str(t_eqSol)]);
%     display(['Dynamics, Computation time (sec): ', num2str(t_dynSol)]);
%     
%     PlotDDFT(v2struct(optsPhys,optsNum,data));    
             
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function dxdt = dx_dt(t,rho_s)
        
        %x       = x(2:end);        
        I       = eye(N1*N2);

        mu_s     = GetExcessChemPotential(rho_s,t,mu); 
        %h_s      = Diff.grad*rho_s - Vext_grad;
                
        dxdt     = kBT*Diff.Lap*mu_s;% + [I I]*(h_s.*(Diff.grad*mu_s));  
        
        %Boundary Conditions: no flux at the walls        
        flux_dir        = Diff.grad*mu_s;
        dxdt(Ind.bound) = Ind.normal*flux_dir;           
        dxdt(Ind.left)  = Diff.Dy2(Ind.left,:)*rho_s; %BC at origin
             
        %dxdt = [Int_of_path*GetFlux(x,t) ;dxdt];
    end
    function y = f_eq(rho_s)
        %solves for T*log*rho + Vext                        
        y            = GetExcessChemPotential(rho_s,0,mu); 
        
        %Boundary Conditions
        rhoB         = rho_s(Ind.bound);
        y(Ind.bound) = Cn*(Ind.normal*(Diff.grad*rho_s)) + cos(theta)*(rhoB-1).*rhoB;
        y(Ind.right) = Ind.normalRight*(Diff.grad*rho_s);
        y(Ind.left)  = Diff.Dy2(Ind.left,:)*rho_s; %BC at origin
        %Here, we require continuity of the first derivative. Is this cheating??       
        y(Ind.left&Ind.bottom) = Diff.Dy1(Ind.left&Ind.bottom,:)*rho_s + Diff.Dy1(Ind.left&Ind.top,:)*rho_s;        
        
        %y            = [Int*rho_s - nParticles;y];
    end
    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = GetExcessChemPotential(x,t,mu); 
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end

    function mu_s = GetExcessChemPotential(rho_s,t,mu_offset)    
        [fE,dfE]    = f(rho_s);
        mu_s        = (rho_s.*dfE + fE)/Cn - Cn*(Diff.Lap*rho_s) + ...
                        - mu_offset + ...
                        getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys);                
    end   
    function [z,dz] = f(rho)
        z   = 0.5*rho.*(1-rho).^2;
        dz  = 0.5*(1-rho).*(1-3*rho);
    end
    function z = InitialGuess()
        xs = Pts.y1_kv.*cos(Pts.y2_kv);
        ys = Pts.y1_kv.*sin(Pts.y2_kv);
        z = 0.5*(1 - tanh((xs*cos(theta)-ys*sin(theta))/(2*Cn)));
    end


end