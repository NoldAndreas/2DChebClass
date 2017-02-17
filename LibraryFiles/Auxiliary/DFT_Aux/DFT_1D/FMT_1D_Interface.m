function [rho_ic1D,postParms] = FMT_1D_Interface(HS,IntMatrFex_2D,optsPhys,FexNum,Conv,BoolPlot,yShift)

    saveFigs = true;
    
    global dirData
    
    %************************************************
    %***************  Initialization ****************
    %************************************************
    kBT       = optsPhys.kBT;           
    R         = optsPhys.sigmaS/2;
    fBulk     = str2func(['FexBulk_',FexNum.Fex]);  
    
    if(nargin <7)
        yShift = 0;
    end    
    if(nargin < 6)
        BoolPlot = true;
    end
    
    if(isfield(optsPhys,'Dmu') && isfield(optsPhys,'mu_sat'))
        mu        = optsPhys.mu_sat + optsPhys.Dmu;
    elseif(isfield(optsPhys,'eta'))
        eta       = optsPhys.eta;
        rhoBulk   = eta*6/pi;
        mu        = kBT*log(rhoBulk) + fBulk(rhoBulk,kBT);
    end
    
    getFex = str2func(['Fex_',FexNum.Fex]);
    
    markComp  = (HS.Pts.y2_kv==inf);     
    nSpecies  = 1;
    Pts       = HS.Pts;
    
    %************************************************
    %****************  Preprocess  ******************
    %************************************************        
    markADComp     = (HS.AD.Pts.y2_kv==inf); 
    
    
    [h1,Int_1D,h2]     = HS.ComputeIntegrationVector();
        
   if(isempty(IntMatrFex_2D)) %LDA Approaches
        IntMatrFex = [];   
        Int_1D_AD  = [];
    else %FMT approaches:
        [h1,Int_1D,h2]     = HS.ComputeIntegrationVector();
        [h1,Int_1D_AD,h2]  = HS.AD.ComputeIntegrationVector();
        %[h1,Int_1D_AD,h2]  = HS.AD.ComputeIntegrationVector();        
        %[h1,h2,Int_1D_AD]  = HS.AD.ComputeIntegrationVector();
        
        IntMatrFex_1D  = Get1DMatrices(IntMatrFex_2D,HS,markComp,markADComp);
        IntMatrFex     = IntMatrFex_1D;
    end

     if(isempty(Conv))
        Conv = zeros(size(IntMatrFex));
    else
        Conv = Conv(markComp,markComp);
    end    
           
    y1S     = Pts.y1_kv(markComp); 
    y2S     = Pts.y2_kv(markComp);
    	
    if(isa(HS,'HalfSpaceSkewed'))        
        %consider that integration is done normal to interface 
        Int_1D    = Int_1D*sin(HS.alpha);
        Int_1D_AD = Int_1D_AD*sin(HS.alpha);
    end
    
    %****************************************************************
    %**************** Solve for equilibrium 1D condition   **********
    %****************************************************************        
    tic
    if(~isfield(optsPhys,'rho_iguess'))
        optsPhys.rho_iguess = (optsPhys.rhoLiq_sat+optsPhys.rhoGas_sat)/2 + ...
                              (optsPhys.rhoLiq_sat-optsPhys.rhoGas_sat)/2*tanh(y1S-yShift);
        y0 = kBT*log(optsPhys.rho_iguess);
    elseif(length(optsPhys.rho_iguess)==1)
        y0 = kBT*log(optsPhys.rho_iguess)*ones(length(Pts.y1_kv),nSpecies);    
        y0 = y0(markComp);       
    else
        y0 = kBT*log(optsPhys.rho_iguess);                    
    end     
    
    IP0 = barychebevalMatrix(Pts.x1,0);
   
    %PlotRosenfeldFMT_AverageDensities(HS,IntMatrFex(1),ones(size(y0)));                       
    fsolveOpts    = optimset('MaxFunEvals',2000000,'MaxIter',200000,'Algorithm','Levenberg-Marquardt','Display','off');%'TolFun',1e-10,'TolX',1e-10);
    x_ic_1D       = fsolve(@f,y0,fsolveOpts);
    rho_ic1D      = exp(x_ic_1D/kBT);     
    postParms.Fex = GetExcessGrandPotential(rho_ic1D);    

    %*****************************************
    %**************** PostProcess   **********
    %*****************************************
    if(BoolPlot)
        
        f1 = figure('Color','white','Position', [0 0 1000 1000],...
                    'name','Liquid-vapour interface');        
        HS.do1DPlotParallel(rho_ic1D); hold on;
        plot([-5 5],[optsPhys.rhoLiq_sat,optsPhys.rhoLiq_sat],'--','linewidth',2);
        plot([-5 5],[optsPhys.rhoGas_sat,optsPhys.rhoGas_sat],'--','linewidth',2);
        hold off;
                        
        
         if(saveFigs && ~isempty(IntMatrFex))
            nStruct = IntMatrFex(1).AD;
            
            f2 = figure('Color','white');
            set(f2, 'Position', [0 0 400 400]);
            HS.do1DPlotParallel(nStruct.n2*rho_ic1D); 
            hh = title('$\rho_2$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
            pbaspect([1 1 1]);            
            xlim([-4 4]);
            set(gca,'Xtick',[0 4]);
            inset2(f1,f2,0.25,[0.2,0.6]);
            close(f2);            
            
            f2 = figure('Color','white');
            set(f2, 'Position', [0 0 400 400]);
            HS.do1DPlotParallel(nStruct.n3*rho_ic1D); 
            hh = title('$\rho_3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
            pbaspect([1 1 1]);            
            xlim([-4 4]);
            set(gca,'Xtick',[0 4]);
            inset2(f1,f2,0.25,[0.45,0.6]);
            close(f2);            
            
            f2 = figure('Color','white');
            set(f2, 'Position', [0 0 400 400]);
            HS.do1DPlotParallel(nStruct.n2_v_2*rho_ic1D); 
            hh = title('$\rho_3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
            pbaspect([1 1 1]);            
            xlim([-4 4]);
            set(gca,'Xtick',[0 4]);
            inset2(f1,f2,0.25,[0.7,0.6]);
            close(f2);            
            
            colormap(gray);    
            set(gcf,'Color','white');            
            set(gcf, 'Position', [0 0 1400 1000]);
            
            if(~exist(dirData,'dir'))                        
                mkdir(dirData);
            end
            %SaveFigure('Density_Interface');                        
         end
    end
    
    
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************                      
    function y = f(x)
        %solves for T*log*rho + Vext                        
        y            = GetExcessChemPotential(x,0,mu);         
        y            = [y(:).*exp(-x/kBT);...
                        IP0*(x - y0)];
                        %x(ceil(end/2))-y0(ceil(end/2))];
    end
    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp(x/kBT);                
        mu_s  = getFex(rho_s,IntMatrFex,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end
        Vadd = getVAdd(y1S,y2S,t,optsPhys.V1);
        mu_s = mu_s + x + Conv*rho_s + Vadd;
    end


    function Fex = GetExcessGrandPotential(rho)
        %Compute excess grand potential
        f_id       = kBT*rho.*(log(rho)-1);
        [h1s,h2s,f_hs] = getFex(rho,IntMatrFex,kBT,R);
        f_attr     = 0.5*rho.*(Conv*rho);
        f_Vmu      = - mu*rho;        
        
        f_loc      = f_id + f_attr + f_Vmu;      
        
        if(isempty(IntMatrFex))
            Fex    = Int_1D*( (f_hs-f_hs(end)) + (f_loc-f_loc(end)) );
        else
            Fex    = Int_1D_AD*(f_hs-f_hs(end)) + Int_1D*(f_loc-f_loc(end));
        end                
        
        %Get Bulk value and compare:
        floc_Bulk = GetfBulk(rho(end));        
        PrintErrorPos(floc_Bulk-(f_loc(end)+f_hs(end)),'Bulk Free Energy in 1D Interface Computation at inf ');
        if(length(y1S) > 3)
            PrintErrorPos(floc_Bulk-(f_loc(end-1)+f_hs(end-1)),['Bulk Free Energy in 1D Interface Computation at ',num2str(y1S(end-1))]);
        end
        floc_Bulk = GetfBulk(rho(1));
        PrintErrorPos(floc_Bulk-(f_loc(1)+f_hs(1)),'Bulk Free Energy in 1D Interface Computation at -inf ');        
        if(length(y1S) > 3)
            PrintErrorPos(floc_Bulk-(f_loc(2)+f_hs(2)),['Bulk Free Energy in 1D Interface Computation at ',num2str(y1S(2))]);
        end
    end

    function floc_Bulk = GetfBulk(rho_Bulk)        
        [h1s,f_hs_Bulk,h2s,h3s] = FexBulk_FMTRosenfeld_3DFluid(rho_Bulk,kBT);
        Phi_r             = str2func(optsPhys.V2.V2DV2);        
        [h1s,h2s,alpha]   = Phi_r(0,optsPhys.V2);
        f_attr_Bulk       = alpha*rho_Bulk^2;
        f_id_Bulk         = kBT*rho_Bulk.*(log(rho_Bulk)-1);
        f_Vmu_Bulk        = -mu*rho_Bulk;
        
        floc_Bulk         = f_id_Bulk + f_attr_Bulk + f_hs_Bulk + f_Vmu_Bulk;
    end

    %***************************************************************
    %   Plotting Auxiliary functions:
    %***************************************************************     
	function PlotRosenfeldFMT_AverageDensitiesInf(FMTMatrices,rho)
        
        f1 = figure('name','Average densities');
        set(gcf,'Color','white');
        set(f1, 'Position', [0 0 1300 400]);
        
        nStruct = FMTMatrices.AD;            
        
        subplot(1,3,1);            
        do1Dplot(nStruct.n2*rho); 
        hh = title('$n_2$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
        pbaspect([1 1 1]);

        subplot(1,3,2);
        do1Dplot(nStruct.n3*rho);  
        hh = title('$n_3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
        pbaspect([1 1 1]);

        subplot(1,3,3);
        do1Dplot(nStruct.n2_v_2*rho);  
        hh = title('${\bf n}_{2,y}$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);  
        pbaspect([1 1 1]);

    end

    
   function PlotRosenfeldFMT_AADInf(FMTMatrices,rho)
        
        f1 = figure('name','AAD');
        set(gcf,'Color','white');
        set(f1, 'Position', [0 0 1300 400]);
        nStruct = FMTMatrices.AAD;                
        
        subplot(1,3,1); 
        do1Dplot_D(nStruct.n2*rho); 
        hh = title('$n_2$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
        pbaspect([1 1 1]);
        
        subplot(1,3,2); 
        do1Dplot_D(nStruct.n3*rho); title('n3');    
        hh = title('$n_3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
        pbaspect([1 1 1]);
        
        subplot(1,3,3); 
        do1Dplot_D(nStruct.n2_v_2*rho); title('n2_2');    
        hh = title('${\bf n}_{2,y}$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
        pbaspect([1 1 1]);

    end
    function do1Dplot(val)
        HS.AD.do1DPlotParallel(val);
    end

    function do1Dplot_D(val)
        HS.do1DPlotParallel(val);
    end



end