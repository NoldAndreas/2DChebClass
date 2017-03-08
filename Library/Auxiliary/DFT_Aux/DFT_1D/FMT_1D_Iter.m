function [rho_ic1D,postParms] = FMT_1D_Iter(HS,IntMatrFex_2D,optsPhys,FexNum,Conv,opts)
                
    if(nargin < 6)
        opts = {};
    else
        if(islogical(opts))
            if(opts)
                opts = {'plot'};
            else
                opts = {};
            end            
        end
    end
    if(IsOption(opts,'NumericsManuscript'))
        xLabelTxt = '$y_2$'; yLabelTxt = '$n$';
    else
        xLabelTxt = '$y_2/\sigma$'; yLabelTxt = '$n\sigma^3$';
    end        

    %************************************************
    %***************  Initialization ****************
    %************************************************
    kBT             = optsPhys.kBT;           
    R               = optsPhys.sigmaS/2;
    fBulk           = str2func(['FexBulk_',FexNum.Fex]);  
    optsPhys.HSBulk = (['FexBulk_',FexNum.Fex]);
    getFex          = str2func(['Fex_',FexNum.Fex]);
	markComp        = (HS.Pts.y1_kv==inf);        
    nSpecies        = 1;
    Pts             = HS.Pts;
    N               = HS.M;
    N_AD            = HS.AD.M;
    
    if(~isfield(optsPhys,'V2'))
        optsPhys.V2 = struct('V2DV2','zeroPotential');                  
    end
    
    if(isfield(optsPhys,'Dmu') && isfield(optsPhys,'mu_sat'))
        mu        = optsPhys.mu_sat + optsPhys.Dmu;        
    elseif(isfield(optsPhys,'eta'))
        eta       = optsPhys.eta;
        rhoBulk   = eta*6/pi;
        mu        = kBT*log(rhoBulk) + fBulk(rhoBulk,kBT);               
    end    
    
    [rhoGas_eq,rhoLiq_eq,pLiq,pGas] = BulkValues(mu,optsPhys,[],false);    
    if(abs(rhoGas_eq - optsPhys.rho_iguess(end)) < abs(rhoLiq_eq - optsPhys.rho_iguess(end)))
        rhoBulk = rhoGas_eq;
        pBulk   = pGas;
    else
        rhoBulk = rhoLiq_eq;
        pBulk   = pLiq;
    end
                
    %************************************************
    %****************  Preprocess  ******************
    %************************************************	
    y1S     = Pts.y1_kv(markComp); 
    y2S     = Pts.y2_kv(markComp);
    ptsCart = HS.GetCartPts(y1S,y2S);
    
    [VAdd,dVAdd]    = getVAdd(ptsCart.y1_kv,ptsCart.y2_kv,0,optsPhys.V1);    
	[h1,h2,Int_1D]  = HS.ComputeIntegrationVector();    
    
    if(isempty(IntMatrFex_2D)) %LDA Approaches
        IntMatrFex = [];
        Int_1D_AD  = Int_1D;
    else %FMT approaches:
        [h1,h2,Int_1D_AD]  = HS.AD.ComputeIntegrationVector();
        
        IntMatrFex_1D     = Get1DMatrices(IntMatrFex_2D,HS);    
        IntMatrFex        = IntMatrFex_1D;   
    end
    
    if(isempty(Conv))
        Conv = zeros(size(IntMatrFex));
    elseif(length(Conv)~= sum(markComp))
        Conv = Conv(markComp,markComp);
    end                   
    
    y2MaxInt = inf;
	Int_1D(ptsCart.y2_kv>y2MaxInt) = 0;
    mark      = (HS.AD.Pts.y1_kv == inf);
    PtsADCart = HS.AD.GetCartPts();
    Int_1D_AD(PtsADCart.y2_kv(mark)>y2MaxInt) = 0;
    
    if(y2MaxInt < Inf)
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
        y0        = zeros(length(Pts.y1_kv),nSpecies);
    end       
    y0 = y0(markComp);        
    
    if(IsOption(opts,'Newton'))
        [x_ic_1D] = NewtonMethod(zeros(N+4*N_AD,1),@f,1,3,0.3,{'returnLastIteration'});
        [x_ic_1D] = NewtonMethod(x_ic_1D,@f,1e-10,200,1);   
    elseif(IsOption(opts,'Picard'))
        xmin = -10;
        xmax = 20;
        x_ic = y0;
        err  = 1; it = 1;
        while((err > 1e-10) && (it < 4000))
            [dx,err]   = GetdX_nP1(x_ic);
            x_ic       = x_ic + 0.05*dx;
            x_ic       = max(min(x_ic,xmax),xmin);
            disp(['Iteration ',num2str(it),': ',num2str(err)]);
            it = it+1;
        end
        x_ic_1D = x_ic;
	elseif(IsOption(opts,'PicardBound'))
        xmin = kBT*log(0.5*rhoGas_eq);
        xmax = kBT*log(10*rhoLiq_eq);
        x_ic = y0;
        err  = 1; it = 1; 
        while((err > 1e-10) && (it < 600))
                        
            [dx,err]   = GetdX_nP1(x_ic);
            l_max  = (xmax - x_ic)./dx; l_max(l_max < 0) = 100;
            l_min  = (xmin - x_ic)./dx; l_min(l_min < 0) = 100;

            lambda = min(1,0.1*min(l_max,l_min));
            x_ic       = x_ic + lambda.*dx;            
            x_ic(HS.Pts.y2 > 5) = kBT*log(rhoGas_eq);
            disp(['Iteration ',num2str(it),': ',num2str(err),' min(lambda) = ',num2str(min(lambda)),' max(lambda) = ',num2str(max(lambda))]);
            it = it+1;
        end
        x_ic_1D = x_ic;	
    elseif(IsOption(opts,'PicardRho'))  
        rho_min = (0.5*rhoGas_eq);
        rho_max = (7*rhoLiq_eq);
        rho_ic = rhoGas_eq*ones(N,1);
        err  = 1; it = 1; 
        while((err > 1e-10) && (it < 600))                        
            [drho,err]   = GetdRho_nP1(rho_ic);
            drho(HS.Pts.y2 > 10) = 0;
            l_max  = (rho_max - rho_ic)./drho; l_max(l_max < 0) = 100;
            l_min  = (rho_min - rho_ic)./drho; l_min(l_min < 0) = 100;
            lambda_maxAbs = 0.05;
            lambda = min(lambda_maxAbs,0.8*min(l_max,l_min));            
            rho_ic       = rho_ic + lambda.*drho;            
            rho_ic(HS.Pts.y2 > 10) = rhoGas_eq;
            [err,i_err] = max(abs(drho));
            %x_ic       = max(min(x_ic,xmax),xmin);
            disp(['Iteration ',num2str(it),' , error  : ',num2str(err),' at y_2 = ',num2str(HS.Pts.y2(i_err)),' min(lambda) = ',num2str(min(lambda)),' max(lambda) = ',num2str(max(lambda))]);
            it = it+1;
        end
        x_ic_1D = kBT*log(rho_ic);	        
	else
        disp('solver not working in this file yet. Check FMT_1D.m');
        fsolveOpts=optimset('Display','off','TolFun',1e-10,'TolX',1e-10);%'MaxFunEvals',2000000,'MaxIter',200000,...
        [x_ic_1D,h1,flag] = fsolve(@f,y0,fsolveOpts);     
        if(flag ~= 1)
            cprintf('red','Error in fsolve, FMT_1D_Interface');
        end
    end
    
    
    rho_ic1D  = exp(x_ic_1D(1:N)/kBT); 
    
	%*****************************************
    %**************** PostProcess   **********
    %*****************************************
    %Check Contact Density: see also Eq. (13a) of [Swol,Henderson,PRA,Vol 40,2567]    
    if(~isempty(dVAdd.dy2))
        checkContactDensity               = (pBulk + Int_1D*(rho_ic1D.*dVAdd.dy2) )/kBT;
        postParms.contactDensity_relError = ((rho_ic1D(1)-checkContactDensity)/checkContactDensity);
        PrintErrorPos(rho_ic1D(1)-checkContactDensity,'First Sum Rule - for Contact Density');
        PrintErrorPos(postParms.contactDensity_relError*100,'First Sum Rule - for contact density (in per cent)');
    else
        checkContactDensity = [];
    end
    
    %****************************
    %********** Plot ************
    %****************************
    if(IsOption(opts,'plot'))        
        if(IsOption(opts,'NoCollPts'))
            bool_collPts = [];
        else
            bool_collPts = 'o';
        end
        f1 = figure;        
        
        subplot(3,3,[1,2,4,5,7,8]);
        HS.do1DPlotNormal(rho_ic1D,bool_collPts); hold on;
        h = xlabel(xLabelTxt);  set(h,'Interpreter','Latex'); set(h,'fontsize',25);
        h = ylabel(yLabelTxt);  set(h,'Interpreter','Latex'); set(h,'fontsize',25);
        
        plot([0 10],[rhoBulk,rhoBulk],'k--');
        if(~isempty(IntMatrFex) && ~isempty(checkContactDensity))
            plot([0 10],[checkContactDensity,checkContactDensity],'k--');
        end               
        xlim([0 7]);
            
        deltaY = 0.5;
        if(strcmp(optsPhys.V2.V2DV2,'zeroPotential') && strcmp(optsPhys.HSBulk,'FexBulk_FMTRosenfeld_3DFluid'))                
            Interp1D              = HS.ComputeInterpolationMatrix(1,(-1:0.01:0.7)',true,true);        
            Interp1D.InterPol     = Interp1D.InterPol(:,markComp);
            Interp1D.ptsCart      = HS.GetCartPts(Interp1D.pts1,Interp1D.pts2);
            dataMC                = LoadGrootData(eta*6/pi);                        
            plot(dataMC.y-deltaY,dataMC.rho,'ks','markerFace','k'); hold on;
            
            f2 = figure('Color','white');
            plot(dataMC.y-deltaY,dataMC.rho,'ks','markerFace','k'); hold on;
            if(eta == 0.4783)
                xlim([0.5 .6]-0.5);   ylim([0 12]);
            elseif(abs(eta - 0.3744) < 0.0001)
                xlim([0.8 1.4]);   ylim([0.5 1.2]);
            elseif(eta == 0.4257)
                xlim([1.3 1.8]-0.5);   ylim([0.65 1.55]);
            end
            plot(Interp1D.ptsCart.y2_kv-deltaY,Interp1D.InterPol*rho_ic1D,'k');
            h = xlabel(xLabelTxt);  set(h,'Interpreter','Latex'); set(h,'fontsize',25);
            h = ylabel(yLabelTxt);  set(h,'Interpreter','Latex'); set(h,'fontsize',25);
            pbaspect([1 1 1]);                            
            if(~IsOption(opts,'plotTex'))
                set(gca,'fontsize',20); 
                set(gca,'linewidth',1.5);                
            end
            hold off;           

            inset(f1,f2,0.35,[.35 .6]);   colormap(gray);    
            set(gcf,'Color','white');  close(f2); close(f1);
            
            if(IsOption(opts,'plotTex'))
                set(gcf, 'Position', [0 0 400 300]);                        
            else
                set(gcf, 'Position', [0 0 1400 1000]);            
            end
        end

               
         %*************** Plot Average Densities ***************        
        if(~isempty(IntMatrFex)) %PlotRosenfeldFMT_AverageDensitiesInf(IntMatrFex(1),rho_ic1D);                
            
            
            nStruct = IntMatrFex(1).AD; 
            
            %f2 = figure('Color','white');
            %set(f2, 'Position', [0 0 400 400]);
            subplot(3,3,3);
            do1Dplot(nStruct.n2*rho_ic1D,bool_collPts);             
            if(IsOption(opts,'NumericsManuscript'))            
                ylabel('$n_2$','Interpreter','Latex','fontsize',25);
            else
                ylabel('$n_2 \sigma^3$','Interpreter','Latex','fontsize',25);
            end
            xlabel(xLabelTxt,'Interpreter','Latex','fontsize',25);
            pbaspect([1 1 1]);            
            xlim([-0.5 3.5]); set(gca,'Xtick',[-0.5 3.5]);
            %inset2(f1,f2,0.25,[0.2,0.6]); close(f2);            
            
            %f2 = figure('Color','white'); set(f2, 'Position', [0 0 400 400]);
            subplot(3,3,6);
            do1Dplot(nStruct.n3*rho_ic1D,bool_collPts);             
            if(IsOption(opts,'NumericsManuscript'))            
                ylabel('$n_3$','Interpreter','Latex','fontsize',25);
            else
                ylabel('$n_3 \sigma^3$','Interpreter','Latex','fontsize',25);
            end
            xlabel(xLabelTxt,'Interpreter','Latex','fontsize',25);
            pbaspect([1 1 1]);            
            xlim([-0.5 3.5]); set(gca,'Xtick',[-0.5 3.5]);
            %inset2(f1,f2,0.25,[0.45,0.6]); close(f2); 
            
            %f2 = figure('Color','white'); set(f2, 'Position', [0 0 400 400]);
            subplot(3,3,9);
            do1Dplot(nStruct.n2_v_2*rho_ic1D,bool_collPts); 
            if(IsOption(opts,'NumericsManuscript'))            
                ylabel('${\left(\bf n_{2}\right)}_2$','Interpreter','Latex','fontsize',25);
            else
                ylabel('${\left(\bf n_{2}\right)}_2 \sigma^3$','Interpreter','Latex','fontsize',25);
            end            
            %ylabel('${\left(\bf n_{2}\right)}_2 \sigma^3$','Interpreter','Latex','fontsize',25);
            xlabel(xLabelTxt,'Interpreter','Latex','fontsize',25);
            pbaspect([1 1 1]);            
            xlim([-0.5 3.5]); set(gca,'Xtick',[-0.5 3.5]);
            %inset2(f1,f2,0.25,[0.7,0.6]); close(f2);            
            
            colormap(gray);    
            set(gcf,'Color','white');            
        end
        SaveFigure(['Density_Wall_rho=0_',num2str(ceil(100*rho_ic1D(end)))]);                        
        %***************************************************************
        %figure('name','Variation of FMT Excess Free Energy for initial condition');
        %do1Dplot_D(Fex_FMTRosenfeld_3DFluid(rho_ic1D,IntMatrFex,kBT));
        %***************************************************************

    end
    
    
    %***************************************************************
    %   Physical Auxiliary functions:    
    %***************************************************************             
    function [y,J] = f(x)
        %solves for T*log*rho + Vext                        
        [y,J]        = GetExcessChemPotential(x,0,mu);         
        y            = y(:);        
    end
    function [mu_s,J_s] = GetExcessChemPotential(x,t,mu_offset)
        xR              = (x(1:N));        
        rho_s           = exp(xR/kBT);
        x(1:N)          = rho_s;
        [mu_s,~,~,J_s]  = getFex(x,IntMatrFex,kBT,R);        
        
        mu_s(1:N)     = mu_s(1:N) + Conv*rho_s + xR + VAdd - mu_offset(1);   
        
        if(~isempty(J_s))            
            J_s(1:N,1:N)  = (J_s(1:N,1:N) + Conv)*diag(rho_s(:))/kBT + eye(N*nSpecies);      
        end

        %for iSpecies=1:nSpecies
%           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
%        end        
 %       mu_s = mu_s + x + Conv*rho_s + VAdd;
    end
    function Fex = GetExcessGrandPotential(rho)
        %Compute excess grand potential
        f_id       = kBT*rho.*(log(rho)-1);
        [h1s,h2s,f_hs] = getFex(rho,IntMatrFex,kBT,R);
        f_attr         = 0.5*rho.*(Conv*rho);
        f_Vmu          = rho.*(VAdd- mu);        
        
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
        [h1s,h2s,alpha]   = getV2(0,optsPhys.V2);    
        f_attr_Bulk       = alpha*rho_Bulk^2;
        f_id_Bulk         = kBT*rho_Bulk.*(log(rho_Bulk)-1);
        f_Vmu_Bulk        = -mu*rho_Bulk;
        
        floc_Bulk         = f_id_Bulk + f_attr_Bulk + f_hs_Bulk + f_Vmu_Bulk;
        PrintErrorPos(floc_Bulk-(f_loc(end)+f_hs(end)),['Bulk Free Energy in 1D Computation for rho=',num2str(rho(end),2),'at inf']);
        PrintErrorPos(floc_Bulk-(f_loc(end-1)+f_hs(end-1)),...
            ['Bulk Free Energy in 1D Computation for rho=',num2str(rho(end),2),'at ',num2str(y2S(end-1),2)]);
    end

    function [dx,err] = GetdX_nP1(x)
        %x       = GetFullX(x);
        dx      = -GetExcessChemPotential(x,0,mu);    
        %dx      = dx(markFull);
        err     = max(abs(dx));
    end

    function [drho,err] = GetdRho_nP1(rho)
        %x       = GetFullX(x);
        x       = kBT*log(rho);
        Dx      = -GetExcessChemPotential(x,0,mu);    
        rho_n1  = exp((x+Dx)/kBT);        
        drho    = rho_n1-rho;
        %dx      = dx(markFull);
        err     = max(abs(drho));
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
        
        figure('name','AAD','Color','white','Position', [0 0 1300 400]);
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
   function do1Dplot(val,bool_collPts)
       if(nargin == 1)
           bool_collPts = 'o';
       end       
       HS.AD.do1DPlotNormal(val,bool_collPts,-0.5);
   end
   function do1Dplot_D(val,bool_collPts)
       if(nargin == 1)
           bool_collPts = 'o';
       end
        HS.do1DPlotNormal(val,bool_collPts);
   end

    function data = LoadGrootData(rhoB)
        filename = ['Fex' filesep 'FMT' filesep 'SpecialPlotting' filesep 'GrootEtAlData.txt'];                
        
        fid = fopen(filename);
        y = textscan(fid,'',1,'headerlines',3); %[T, rhoG, rhoL]        
        x = textscan(fid,'%f %f %f %f %f %f'); %[T, rhoG, rhoL]
        fclose(fid);         
                
        data = struct('y',[],'rho',[]);
        for i = 2:(length(y))
            if(abs(rhoB - y{i})<1e-4)
                data.rho = x{i};
                data.y   = x{1};
                break;
            end            
        end                
               
    end


end