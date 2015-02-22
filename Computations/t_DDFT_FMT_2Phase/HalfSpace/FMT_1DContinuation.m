 function [rho_cont,ell,par,OmEx,dmuCheck,pts,params] = FMT_1DContinuation(HS,IntMatrFex_2D,optsPhys,FexNum,Conv,parCont,output)
    %parCont = {'mu'}
    %global dirDataOrg
    if((nargin >= 7) && (ischar(output) && ~strcmp(output,'Movie')))
        %ChangeDirData([dirDataOrg filesep 'deg90']);
        ignoreList = {'thetaDirection','rho_iguess'};
        [res,~,params]  = DataStorage('IterativeContinuationPostProcess',...
                        @PostProcessAdsorptionIsotherm,optsPhys,[],output,ignoreList);                  
                    
        %[res,~,params]  = DataStorage('IterativeContinuationPostProcess',...
        %                @PostProcessAdsorptionIsotherm,optsPhys,[],2,ignoreList);                                      
     %   ChangeDirData([dirDataOrg filesep 'deg',num2str(this.optsNum.PhysArea.alpha_deg,3)]);
                            
        par         = res.dmu;
        ell         = res.ell;
        OmEx        = res.OmEx;
        dmuCheck    = res.dmuCheck; 
        if(isfield(res,'pts'))
            pts     = res.pts;
        else
            pts     = [];
        end
        if(isfield(res,'rho_cont'))
            rho_cont  = res.rho_cont;
        else
            rho_cont    = [];
        end        
        return;
    end


    close all;  
    disp('* FMT 1D HardWall Continuation *');
    %*******************************************
    %*************** Parameters ****************
    %*******************************************
    if(~isfield(optsPhys,'NContIterations'))
        NComps = 20;
    else
        NComps = optsPhys.NContIterations;
    end
    
    %************************************************
    %***************  Initialization ****************
    %************************************************
    kBT       = optsPhys.kBT;    
    R         = optsPhys.sigmaS/2;
    if(strcmp(parCont,'mu'))
        par_start  = optsPhys.mu_sat + optsPhys.Dmu;
    elseif(strcmp(parCont,'epw'))
        mu        = optsPhys.mu_sat + optsPhys.Dmu;
        par_start = optsPhys.V1.epsilon_w;
    end
        
    getFex        = str2func(['Fex_',FexNum.Fex]);    
    markComp      = (HS.Pts.y1_kv==inf);   
    
    Pts           = HS.Pts;
    nSpecies      = 1;
    
    y1S     = Pts.y1_kv(markComp); 
    y2S     = Pts.y2_kv(markComp);
    ptsCart = HS.GetCartPts(y1S,y2S);    
    [VAdd,dVAdd] = getVAdd(ptsCart.y1_kv,ptsCart.y2_kv,0,optsPhys.V1);    
    [rhoGas_sat,rhoLiq_sat,mu_sat,p] = BulkSatValues(optsPhys);
    
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
                   
    if(~isempty(HS.Interp))
        y2Max = max(HS.Interp.pts2);
    else
        y2Max = 10;
    end
 	subPts.y2_kv          = (0.5:0.01:y2Max)';
    subPts.y1_kv          = inf*ones(size(subPts.y2_kv));           
    IP                    = HS.SubShapePts(subPts);
    Interp1D.InterPol     = IP(:,markComp);
    Interp1D.pts1         = subPts.y1_kv; 
    Interp1D.pts2         = subPts.y2_kv;
    
	[h1,h2,Int_1D]     = HS.ComputeIntegrationVector();
    [h1,h2,Int_1D_AD]  = HS.AD.ComputeIntegrationVector();
    
    
    opts_h                   = optsPhys;        
	if(strcmp(optsPhys.drying,'drying'))        
        opts_h.rho_iguess        = 'WL';
        optsPhys.thetaDirection  = -1;
    else        
        opts_h.rho_iguess        = 'WG';                        
        optsPhys.thetaDirection  = 1;
    end
    optsPhys.rho_iguess      = FMT_1D(HS,IntMatrFex_2D,opts_h,FexNum,Conv,false);    
   
    
    %****************************************************************
    %**************** Solve for equilibrium 1D condition   **********
    %****************************************************************        
    tic
    if(isfield(optsPhys,'rho_iguess'))
        
        y0 = kBT*log(optsPhys.rho_iguess);
        if(length(y0) ~= sum(markComp))
            y0 = y0(markComp);
        end
    else
        y0        = zeros(sum(markComp),1);%getInitialGuess(VAdd0);                
    end       
    
       
    fprintf(1,'Computing continuation ...');    
       
    if(strcmp(parCont,'mu'))
         x_cont  = IterativeContinuation(@f,NComps,optsPhys.thetaDirection*0.01,[par_start;y0],optsPhys,@US_ScalarProduct,@X_to_rho);
    elseif(strcmp(parCont,'epw'))
        x_cont  = IterativeContinuation(@f_epw,NComps,optsPhys.thetaDirection*0.01,[par_start;y0],optsPhys,@US_ScalarProduct);
    end    
    
    %x_cont = BasicContinutation(@f,mu,y0);
    
    %************************************
    %*********** Postprocess   **********
    %************************************        
    Int_1D(y2S > 45) = 0;    
    [res,~,params]  = DataStorage('IterativeContinuationPostProcess',...
                        @PostProcessAdsorptionIsotherm,optsPhys,x_cont);      
                    
    rho_cont    = exp(x_cont(:,2:end)/kBT); 	
    par         = res.dmu;
	ell         = res.ell;
	OmEx        = res.OmEx;
	dmuCheck    = res.dmuCheck;    
    pts         = HS.Pts;
        
    %Check for adsorption isotherm:
    %Check Contact Density: see also Eq. (13a) of [Swol,Henderson,PRA,Vol 40,2567]
%     for i=1:no_cont
%         rho_i               = rho_cont(i,:)';    
%     	[h1,h2,p]           = BulkValues(mu,optsPhys,[],false);    
%         checkContactDensity = (p + Int_1D*(rho_i.*dVAdd.dy2) )/kBT;        
%         PrintErrorPos(rho_ic1D(1)-checkContactDensity,'First Sum Rule - for Contact Density');
%     end
    
    %************************************
    %*********** Plot   **********
    %************************************
   % no_cont = 37;
    gifFile = getMovieFile('Movie_FMT_1DContinuation');           
    
   % par      = par(1:no_cont-1);
%    ell      = ell(1:no_cont-1);
%    OmEx     = OmEx(1:no_cont-1);
%    dmuCheck = dmuCheck(1:no_cont-1);
%    rho_cont = rho_cont(1:no_cont-1,:);
    
     %Check asymptotic behavior for attractive fluid
    if(strcmp(optsPhys.V1.V1DV1,'Vext_Cart_8') && strcmp(optsPhys.V2.V2DV2,'Phi2DLongRange'))
        mu = par;
        
        %Compute Hamaker constant
        DeltaRho  = optsPhys.rhoLiq_sat - optsPhys.rhoGas_sat;
        A = 4*pi^2*DeltaRho*(optsPhys.rhoLiq_sat*optsPhys.V2.epsilon - optsPhys.V1.epsilon_w);
        
        deltaMu = abs(mu - optsPhys.mu_sat);
        ell_Ana = (abs(A)./(6*pi*DeltaRho*deltaMu)).^(1/3);
        
        figure('Color','white','Position',[0 0 1200 600]);         
        loglog(deltaMu,ell,'linewidth',1.5); hold on;                        
        loglog(deltaMu,ell_Ana,'r--','linewidth',1.5);
        
        %Enforce non-scientific ticks for axis
        set(gca,'yticklabel',num2str(get(gca,'xtick')'))
        set(gca,'yticklabel',num2str(get(gca,'ytick')'))
    end
    
    %Movie
    
    if( (nargin >= 7) && strcmp(output,'Movie'))
        figure('Color','white','Position',[0 0 1200 600]);  
        for i=1:length(ell)    
            subplot(2,2,1);
            HS.do1DPlotNormal(rho_cont(i,:));
            %do1Dplot(rho_cont(i,:)');
            ylim([0 max(max(rho_cont))]);        
            ylabel('$\varrho$','Interpreter','Latex','fontsize',25);        
            title(['par = ',num2str(par(i))]);
            drawnow;            
            hold off;

            subplot(2,2,2);
            do1DXplot(rho_cont(i,:)');
            ylim([0 max(max(rho_cont))]);        
            ylabel('$\varrho$','Interpreter','Latex','fontsize',25);        
            title(['par = ',num2str(par(i))]);
            drawnow;            
            hold off;

            subplot(2,2,3);
            plot([optsPhys.mu_sat optsPhys.mu_sat],[min(ell) max(ell)],'r--','linewidth',1.5); hold on;
            plot(par,ell,'linewidth',1.5); hold on;
            plot(par(i),ell(i),'om','MarkerFace','g','MarkerSize',10); hold off;
            title('Adsorption Diagram');
            xlabel('$\mu$','Interpreter','Latex','fontsize',25);
            ylabel('$\ell$','Interpreter','Latex','fontsize',25);
            set(gca,'fontsize',20);
            set(gca,'linewidth',1.5);
            xlim([(min(par)-0.01) (max(par) + 0.01)]);
            ylim([(min(ell)-0.1) (max(ell)+0.1)]);

%             subplot(2,2,4);
%             plot([optsPhys.mu_sat optsPhys.mu_sat],[min(OmEx) max(OmEx)],'r--','linewidth',1.5); hold on;
%             plot(ell,OmEx,'linewidth',1.5); hold on;
%             plot(ell(i),OmEx(i),'om','MarkerFace','g','MarkerSize',10); hold off;
%             title('Adsorption Diagram');
%             xlabel('$\ell$','Interpreter','Latex','fontsize',25);
%             ylabel('$F_{ex}$','Interpreter','Latex','fontsize',25);
%             set(gca,'fontsize',20);
%             set(gca,'linewidth',1.5);
%             xlim([(min(ell)-0.01) (max(ell) + 0.01)]);
%             ylim([(min(OmEx)-0.1) (max(OmEx)+0.1)]);


            Record(i,gifFile,0.1);               
        end    
        disp(['Movie` saved in: ',gifFile]);
    else        
        figure;
        subplot(1,2,1);
        plot([optsPhys.mu_sat optsPhys.mu_sat],[min(ell) max(ell)],'r--','linewidth',1.5); hold on;
        plot(par,ell,'linewidth',1.5); hold on;        
        plot(mu_sat + dmuCheck,ell,'b--','linewidth',1.5); hold on;   
        title('Adsorption Diagram');
        xlabel('$\mu$','Interpreter','Latex','fontsize',25);
        ylabel('$\ell$','Interpreter','Latex','fontsize',25);
        set(gca,'fontsize',20);
        set(gca,'linewidth',1.5);
        xlim([(min(par)-0.01) (max(par) + 0.01)]);
        ylim([(min(ell)-0.1) (max(ell)+0.1)]);

        subplot(1,2,2);
        plot([optsPhys.mu_sat optsPhys.mu_sat],[min(OmEx) max(OmEx)],'r--','linewidth',1.5); hold on;
        plot(ell,OmEx,'linewidth',1.5); hold on;        
        title('Adsorption Diagram');
        xlabel('$\ell$','Interpreter','Latex','fontsize',25);
        ylabel('$F_{ex}$','Interpreter','Latex','fontsize',25);
        set(gca,'fontsize',20);
        set(gca,'linewidth',1.5);
        xlim([(min(ell)-0.01) (max(ell) + 0.01)]);
        ylim([(min(OmEx)-0.1) (max(OmEx)+0.1)]);
    end
    
    %*****************************************
    %**************** PostProcess   **********
    %*****************************************


    %************************************************** *************
    %   Physical Auxiliary functions:
    %***************************************************************             
    function y = f(x,mus)
        %solves for T*log*rho + Vext                        
        y = GetExcessChemPotential(x,0,mus);                
    end

    function y = f_epw(x,epw)
        %solves for T*log*rho + Vext
        optsPhys.V1.epsilon_w = epw;        
        y    = GetExcessChemPotential(x,0,mu);
    end

    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp(x/kBT);                
        mu_s  = getFex(rho_s,IntMatrFex,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end
        Vadd = getVAdd(ptsCart.y1_kv,ptsCart.y2_kv,0,optsPhys.V1);        
        mu_s = mu_s + x + Conv*rho_s + Vadd;
    end

    function OmEx = GetExcessGrandPotential(rho,par)
        %Compute excess grand potential
        f_id       = kBT*rho.*(log(rho)-1);
        [h1s,h2s,f_hs] = getFex(rho,IntMatrFex,kBT,R);
        f_attr     = 0.5*Conv*rho;
        if(strcmp(parCont,'mu'))
            f_Vmu      = rho.*(getVAdd(y1S,y2S,0,optsPhys.V1)- par);
        elseif(strcmp(parCont,'epw'))
            optsPhys.V1.epsilon_w = par;
            f_Vmu      = rho.*(getVAdd(y1S,y2S,0,optsPhys.V1)- mus);
        end
        
        f_loc  = f_id + f_attr + f_Vmu;                 
        OmEx    = Int_1D_AD*(f_hs-f_hs(end)) + Int_1D*(f_loc-f_loc(end));
    end

     function n = US_ScalarProduct(x,y)
         int       = zeros(size(Int_1D));
         mark      = (HS.Pts.y2 < 10);
         int(mark) = Int_1D(mark);
         x = x(:);
         %n = (x(1)).*y(1) + int*((x(2:end)-x(end)).*(y(2:end)-y(end)));
         %n = (x(1)).*y(1) + int*((x(2:end)).*(y(2:end)));
         n = (x(1)).*y(1) + (int*x(2:end)).*(int*y(2:end));
     end
 
     function rhoX = X_to_rho(x)
         rhoX        = x;
         rhoX(2:end) = exp(x(2:end)/kBT);
     end

    %***************************************************************
    %   Plotting Auxiliary functions:
    %***************************************************************     	
    function do1Dplot(val)
            plot(HS.Pts.y2_kv(HS.Pts.y1_kv == inf),val,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); hold on
            plot(Interp1D.pts2,Interp1D.InterPol*val,'linewidth',1.5);
            xlim([min(Interp1D.pts2) max(Interp1D.pts2)]);    
            xlabel('$y_2$','Interpreter','Latex','fontsize',25);                    
        	%h = ylabel('$\rho$');  set(h,'Interpreter','Latex'); set(h,'fontsize',25);                        
            %title(sel(i));    
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);                        
    end

    function do1DXplot(val)
            plot(HS.Pts.x2,val,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); hold on
            plot(HS.CompSpace2(Interp1D.pts2),Interp1D.InterPol*val,'linewidth',1.5);
            xlim([-1 1]);
            xlabel('$x$','Interpreter','Latex','fontsize',25);
        	%h = ylabel('$\rho$');  set(h,'Interpreter','Latex'); set(h,'fontsize',25);                        
            %title(sel(i));    
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);                        
    end


     function res = PostProcessAdsorptionIsotherm(optsPhys,x_cont)        
        rho_cont_h    = exp(x_cont(:,2:end)/kBT); 
        par_h         = x_cont(:,1);    
        
        res.ell       = zeros(size(par_h));
        res.OmEx      = zeros(size(par_h));
        res.dmuCheck  = zeros(size(par_h));                    
        res.dmu       = par_h;
        res.rho_cont  = rho_cont_h;
        res.pts       = HS.Pts;
        
        no_cont_h     = size(x_cont,1);
        
        
        for j=1:no_cont_h
            rho_i      = rho_cont_h(j,:)';
            res.ell(j) = Int_1D/(optsPhys.rhoLiq_sat-optsPhys.rhoGas_sat)*...
                                            (rho_i - rho_i(end));        

            %if( (nargin >= 7) && strcmp(output,'full'))
            res.OmEx(j)    = GetExcessGrandPotential(rho_i,par_h(j));               
            [rhoGas_eq,rhoLiq_eq,pL,pG] = BulkValues(par_h(j),optsPhys,[],false);    
            %Check Sum rule - jn case that par_h = chemical potential
            %(1) Contact Density
            checkContactDensity = (pG + Int_1D*(rho_i.*dVAdd.dy2) )/kBT;
            PrintErrorPos(checkContactDensity-rho_i(1),'contact density');

            %(2) Disjoining Pressure - Following Sum Rule
            % - (mu - mu_sat) Delta(n) = (n(0,h) - n(0,infty)) - int( (n(z,h)  - n(z,infty))*dVext/dz  ,zt = 0 ... infty)
            optss              = optsPhys;
            optss.Dmu          = par_h(j)- mu_sat;
            
            %Decide whether its wetting or drying    
            if(strcmp(optsPhys.drying,'drying'))
                optss.rho_iguess   = 'WG';
                deltaRho           = -(rhoLiq_eq-rhoGas_eq);
            else
                optss.rho_iguess   = 'WL';
                deltaRho           = (rhoLiq_eq-rhoGas_eq);
            end
            
            [rho_l_inf,postParms] = FMT_1D(HS,IntMatrFex_2D,optss,FexNum,Conv,false);

            %optss.rho_iguess  = 'WG';
            %[rho_wg,postParms] = FMT_1D(HS,IntMatrFex_2D,optss,FexNum,Conv,false);                        

            res.dmuCheck(j) = -(kBT*(rho_i(1) - rho_l_inf(1)) - Int_1D*( (rho_i - rho_l_inf).*dVAdd.dy2))/deltaRho;
          %  end
        end   
     end

end