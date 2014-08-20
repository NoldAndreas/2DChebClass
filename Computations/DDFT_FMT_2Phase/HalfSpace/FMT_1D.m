function [rho_ic1D,postParms] = FMT_1D(HS,IntMatrFex_2D,optsPhys,FexNum,Conv,BoolPlot)
    global PersonalUserOutput dirData MinimalOutput
            
    saveFigs = true;
    if(nargin < 6)
        BoolPlot = true;
    end

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
    
    optsPhys.HSBulk = (['FexBulk_',FexNum.Fex]);
    [rhoGas_eq,rhoLiq_eq,p] = BulkValues(mu,optsPhys,[],false);    
    if(abs(rhoGas_eq - optsPhys.rho_iguess(end)) < abs(rhoLiq_eq - optsPhys.rho_iguess(end)))
        rhoBulk = rhoGas_eq;
    else
        rhoBulk = rhoLiq_eq;
    end
    
    getFex = str2func(['Fex_',FexNum.Fex]);
    
    markComp  = (HS.Pts.y1_kv==inf);        
    nSpecies  = 1;
    Pts       = HS.Pts;
    %************************************************
    %****************  Preprocess  ******************
    %************************************************
	
    y1S     = Pts.y1_kv(markComp); 
    y2S     = Pts.y2_kv(markComp);
    ptsCart = HS.GetCartPts(y1S,y2S);
    
    [VAdd,dVAdd] = getVAdd(ptsCart.y1_kv,ptsCart.y2_kv,0,optsPhys.V1);    
    
	[h1,h2,Int_1D]     = HS.ComputeIntegrationVector();    
    
    
    if(isempty(IntMatrFex_2D)) %LDA Approaches
        IntMatrFex = [];
        Int_1D_AD  = Int_1D;
    else %FMT approaches:
        [h1,h2,Int_1D_AD]  = HS.AD.ComputeIntegrationVector();
        
        IntMatrFex_1D  = Get1DMatrices(IntMatrFex_2D,HS);    
        IntMatrFex     = IntMatrFex_1D;
    end
    
    if(isempty(Conv))
        Conv = zeros(size(IntMatrFex));
    elseif(length(Conv)~= sum(markComp))
        Conv = Conv(markComp,markComp);
    end                   
    
    y2MaxInt = 40;
	Int_1D(ptsCart.y2_kv>y2MaxInt) = 0;
    mark      = (HS.AD.Pts.y1_kv == inf);
    PtsADCart = HS.AD.GetCartPts();
    Int_1D_AD(PtsADCart.y2_kv(mark)>y2MaxInt) = 0;
    if(~MinimalOutput)
        cprintf('-k',['For integration, values of rho for y2Cart > ',num2str(y2MaxInt),' are ignored.\n']);
    end
        
    %****************************************************************
    %**************** Solve for equilibrium 1D condition   **********
    %****************************************************************        
    tic
    if(isfield(optsPhys,'rho_iguess'))        
        if(ischar(optsPhys.rho_iguess))
            if(strcmp(optsPhys.rho_iguess,'WG'))
                y0        = kBT*log(rhoGas_eq)*ones(length(Pts.y1_kv),nSpecies);
            elseif(strcmp(optsPhys.rho_iguess,'WL'))
                y0        = kBT*log(rhoLiq_eq)*ones(length(Pts.y1_kv),nSpecies);
            end
        else
            y0        = kBT*log(optsPhys.rho_iguess)*ones(length(Pts.y1_kv),nSpecies);                          
        end
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
    
	%*****************************************
    %**************** PostProcess   **********
    %*****************************************
    postParms.Fex = GetExcessGrandPotential(rho_ic1D);   
    
    %Check Contact Density: see also Eq. (13a) of [Swol,Henderson,PRA,Vol 40,2567]
    checkContactDensity = (p + Int_1D*(rho_ic1D.*dVAdd.dy2) )/kBT;
    %checkContactDensity = (p)/kBT;
    PrintErrorPos(rho_ic1D(1)-checkContactDensity,'First Sum Rule - for Contact Density');
    
    %****************************
    %********** Plot   **********
    %****************************
    if(BoolPlot && (PersonalUserOutput))
        f1 = figure('Color','white');
        set(f1, 'Position', [0 0 1000 1000]);
        HS.do1DPlotNormal(rho_ic1D); hold on;
        
        plot([0 10],[rhoBulk,rhoBulk],'--','linewidth',2);
        if(~isempty(IntMatrFex))
            plot([0 10],[checkContactDensity,checkContactDensity],'--','linewidth',2);
        end
        hold off;

        %*************** Plot Average Densities ***************        
        if(saveFigs && ~isempty(IntMatrFex))
            nStruct = IntMatrFex(1).AD; 
            
            f2 = figure('Color','white');
            set(f2, 'Position', [0 0 400 400]);
            do1Dplot(nStruct.n2*rho_ic1D,false); 
            hh = title('$n_2\sigma^3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
            pbaspect([1 1 1]);            
            xlim([0 4]);
            set(gca,'Xtick',[0 4]);
            inset2(f1,f2,0.25,[0.2,0.6]);
            close(f2);            
            
            f2 = figure('Color','white');
            set(f2, 'Position', [0 0 400 400]);
            do1Dplot(nStruct.n3*rho_ic1D,false); 
            hh = title('$n_3 \sigma^3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
            pbaspect([1 1 1]);            
            xlim([0 4]);
            set(gca,'Xtick',[0 4]);
            inset2(f1,f2,0.25,[0.45,0.6]);
            close(f2);            
            
            f2 = figure('Color','white');
            set(f2, 'Position', [0 0 400 400]);
            do1Dplot(nStruct.n2_v_2*rho_ic1D,false); 
            hh = title('$n_3 \sigma^3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
            pbaspect([1 1 1]);            
            xlim([0 4]);
            set(gca,'Xtick',[0 4]);
            inset2(f1,f2,0.25,[0.7,0.6]);
            close(f2);            
            
            colormap(gray);    
            set(gcf,'Color','white');            
            
            filename = ['Density_Wall_rho=0_',num2str(ceil(100*rho_ic1D(end)))];
            
            if(~exist(dirData,'dir'))                        
                mkdir(dirData);
            end
            
            print2eps([dirData filesep filename],gcf);    
            saveas(gcf,[dirData filesep filename '.fig']);                
        
        elseif(~isempty(IntMatrFex))
            PlotRosenfeldFMT_AverageDensitiesInf(IntMatrFex(1),rho_ic1D);    
        end        
        
        %***************************************************************
        %figure('name','Variation of FMT Excess Free Energy for initial condition');
        %do1Dplot_D(Fex_FMTRosenfeld_3DFluid(rho_ic1D,IntMatrFex,kBT));
        %***************************************************************

    end
    
    
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************             
    function y = f(x)
        %solves for T*log*rho + Vext                        
        y            = GetExcessChemPotential(x,0,mu);         
        y            = y(:);
    end
    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp(x/kBT);                
        mu_s  = getFex(rho_s,IntMatrFex,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end        
        mu_s = mu_s + x + Conv*rho_s + VAdd;
    end
    function Fex = GetExcessGrandPotential(rho)
        %Compute excess grand potential
        f_id       = kBT*rho.*(log(rho)-1);
        [h1s,h2s,f_hs] = getFex(rho,IntMatrFex,kBT,R);
        f_attr     = 0.5*rho.*(Conv*rho);
        f_Vmu      = rho.*(VAdd- mu);        
        
        f_loc  = f_id + f_attr + f_Vmu;   
        
        %Here, we assume that density converges fast enough!        
        %y2MaxInt = 100;
        %Int_1D(ptsCart.y2_kv>y2MaxInt) = 0;
        %PtsADCart = HS.AD.GetCartPts();
        %mark      = (HS.AD.Pts.y1_kv == inf);
        %Int_1D_AD(PtsADCart.y2_kv(mark)>y2MaxInt) = 0;
        if(isempty(IntMatrFex))
            Fex    = Int_1D*( (f_hs-f_hs(end)) + (f_loc-f_loc(end)) );
        else
            Fex    = Int_1D_AD*(f_hs-f_hs(end)) + Int_1D*(f_loc-f_loc(end)) +R*f_hs(end);
        end        
        
        %Get Bulk value and compare:
        rho_Bulk          = rho(end);        
        [h1s,f_hs_Bulk,h2s,h3s] = fBulk(rho_Bulk,kBT);%FexBulk_FMTRosenfeld_3DFluid(rho_Bulk,kBT);
        Phi_r             = str2func(optsPhys.V2.V2DV2);        
        [h1s,h2s,alpha]   = Phi_r(0);    
        f_attr_Bulk       = alpha*rho_Bulk^2;
        f_id_Bulk         = kBT*rho_Bulk.*(log(rho_Bulk)-1);
        f_Vmu_Bulk        = -mu*rho_Bulk;
        
        floc_Bulk         = f_id_Bulk + f_attr_Bulk + f_hs_Bulk + f_Vmu_Bulk;
        PrintErrorPos(floc_Bulk-(f_loc(end)+f_hs(end)),['Bulk Free Energy in 1D Computation for rho=',num2str(rho(end),2),'at inf']);
        PrintErrorPos(floc_Bulk-(f_loc(end-1)+f_hs(end-1)),...
            ['Bulk Free Energy in 1D Computation for rho=',num2str(rho(end),2),'at ',num2str(y2S(end-1),2)]);
    end

    %***************************************************************
    %   Plotting Auxiliary functions:
    %***************************************************************     
	function PlotRosenfeldFMT_AverageDensitiesInf(FMTMatrices,rho)
        
        fs = figure('name','Average densities');
        set(gcf,'Color','white');
        set(fs, 'Position', [0 0 1300 400]);
        
        nStruct = FMTMatrices.AD;            
        
        subplot(1,3,1);            
        do1Dplot(nStruct.n2*rho); 
        hh = title('$n_2 \sigma^3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);        
        xlabel('$y_2/\sigma$','Interpreter','Latex');
        pbaspect([1 1 1]);

        subplot(1,3,2);
        do1Dplot(nStruct.n3*rho);  
        hh = title('$n_3 \sigma^3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);
        xlabel('$y_2/\sigma$','Interpreter','Latex');
        pbaspect([1 1 1]);

        subplot(1,3,3);
        do1Dplot(nStruct.n2_v_2*rho);  
        hh = title('${\bf n}_{2,y} \sigma^3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',25);  
        xlabel('$y_2/\sigma$','Interpreter','Latex');
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
   function do1Dplot(val,BoolSC)
       HS.AD.do1DPlotNormal(val);
   end
   function do1Dplot_D(val)
            HS.do1DPlotNormal(val);
    end



end