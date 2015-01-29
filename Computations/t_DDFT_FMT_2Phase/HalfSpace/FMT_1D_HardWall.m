function rho_ic1D = FMT_1D_HardWall(HS,IntMatrFex_2D,optsPhys,optsNum)
%************************************************************************* 
% define
%   mu_s_i := kBT*log(rho_i) + sum( int(rho_j(r')*Phi2D(r-r'),dr') , 
%               j = 1..N) +mu_HS(rho_i) + V_ext_i - mu_i
% Equilibrium, (solve for each species i)
%  (EQ i,1) mu_s_i     = 0
%  (EQ i,2) int(rho_i) = NoParticles_i
% Dynamics: 
%*************************************************************************       
    global PersonalUserOutput
    %********************************************
    %**************** Initialization   **********
    %********************************************
    
    saveFigs  = true;
    eta       = optsPhys.eta;    
    
    markComp  = (HS.Pts.y1_kv==inf);    
	Pts                       = HS.Pts;
    PtsCart                   = HS.GetCartPts();
    Interp1D                  = HS.ComputeInterpolationMatrix(1,(-1:0.01:0.7)',true,true);        
    Interp1D.InterPol         = Interp1D.InterPol(:,markComp);
    Interp1D.ptsCart          = HS.GetCartPts(Interp1D.pts1,Interp1D.pts2);
         
    subPts.y2_kv              = (0.:0.01:3.5)';
    subPts.y1_kv              = inf*ones(size(subPts.y2_kv));           
    IP                        = HS.AD.SubShapePts(subPts);
    Interp1D_AD.InterPol      = IP(:,HS.AD.Pts.y1_kv == inf);
    Interp1D_AD.pts1          = subPts.y1_kv; 
    Interp1D_AD.pts2          = subPts.y2_kv;
    
	subPtsAAD.y2_kv           = (0.5:0.01:3.5)';
    subPtsAAD.y1_kv           = inf*ones(size(subPtsAAD.y2_kv));           
    IP                        = HS.SubShapePts(subPtsAAD);
    Interp1D_AAD.InterPol     = IP(:,markComp);
    Interp1D_AAD.pts1         = subPtsAAD.y1_kv; 
    Interp1D_AAD.pts2         = subPtsAAD.y2_kv;

	%***********************************
    %**************** Compute **********
    %***********************************
    rho_ic1D = FMT_1D(HS,IntMatrFex_2D,optsPhys,optsNum.FexNum,[],true);
    
	%****************************************************************
    %**************** Postprocess   **********
    %****************************************************************
         
    if(~PersonalUserOutput)
        return;
    end
    %**************** Plot Density ****************
    f1 = figure;
    set(gcf,'Color','white');
    set(f1, 'Position', [0 0 1000 1000]);	
    if(eta == 0.4783)
        shift = struct('xmin',0.235,'xmax',3.31,'ymin',-1.75,'ymax',12.25,'yref',0,'xref',0.5);
        PlotBackgroundImage(['Fex' filesep 'FMT' filesep 'SpecialPlotting' filesep 'PackingFraction0_4783_RothFMTReview.gif'],shift);
        yl = [0 12];
        strEta = '0_4783';
    elseif(eta == 0.4257)
        shift = struct('xmin',0.5,'xmax',3.,'ymin',0.,'ymax',7.,'yref',0,'xref',0.5);
        PlotBackgroundImage(['Fex' filesep 'FMT' filesep 'SpecialPlotting' filesep 'PackingFraction0_4257_RothFMTReview.gif'],shift);
        yl = [0 7];
        strEta = '0_4257';
    end
	%plot(Pts.y2_kv(markComp),rho_ic1D,'o','markersize',1.5); hold on
	%plot(Interp1D.pts2,Interp1D.InterPol*rho_ic,'linewidth',1.5);
	%xlim([0.5 3.5]);                
    
    plot(PtsCart.y2_kv(markComp),rho_ic1D,'o','markersize',8,'markerFace','green'); hold on
    plot(Interp1D.ptsCart.y2_kv,Interp1D.InterPol*rho_ic1D,'linewidth',1.5);
    xlim([0.5 3.]);   ylim(yl);
    h = xlabel('$y/\sigma$');  set(h,'Interpreter','Latex'); set(h,'fontsize',25);
	h = ylabel('$n \sigma^3$');  set(h,'Interpreter','Latex'); set(h,'fontsize',25);
	pbaspect([1 1 1]);                
	set(gca,'fontsize',20);                        
    set(gca,'linewidth',1.5);        
    ax=get(f1,'Position');    
    hold off;
    if(saveFigs)
        print2eps(['Density_eta',strEta],gcf);
        save2pdf(['Density_eta_',num2str(eta),'.pdf'],gcf);
    end

    %**************** Plot Density - Zoom in Close to Wall ****************
	f2 = figure;
    set(gcf,'Color','white');
    if(eta == 0.4783)
        shift = struct('xmin',0.5,'xmax',.6,'ymin',0.,'ymax',11.,'yref',0,'xref',0.5);
        PlotBackgroundImage(['Fex' filesep 'FMT' filesep 'SpecialPlotting' filesep 'Inset_PackingFraction0_4783_RothFMTReview.gif'],shift);
        xlim([0.5 .6]);   ylim([0 12]);
	elseif(eta == 0.4257)
        shift = struct('xmin',1.3,'xmax',1.8,'ymin',0.65,'ymax',1.55,'yref',0.65,'xref',1.3);
        PlotBackgroundImage(['Fex' filesep 'FMT' filesep 'SpecialPlotting' filesep 'Inset_PackingFraction0_4257_RothFMTReview.gif'],shift);
        xlim([1.3 1.8]);   ylim([0.65 1.55]);
    end

	plot(PtsCart.y2_kv(markComp),rho_ic1D,'o','markersize',8,'markerFace','green'); hold on
    plot(Interp1D.ptsCart.y2_kv,Interp1D.InterPol*rho_ic1D,'linewidth',1.5);  %Interp1D.pts2
    h = xlabel('$y/\sigma$');  set(h,'Interpreter','Latex'); set(h,'fontsize',15);
	h = ylabel('$n/\sigma^3$');  set(h,'Interpreter','Latex'); set(h,'fontsize',15);
	pbaspect([1 1 1]);                
	set(gca,'fontsize',15);                        
    set(gca,'linewidth',1.5);                
    hold off;
    
    %save2pdf(['Inset_Density_eta_',num2str(eta),'.pdf'],gcf);   
    
    inset(f1,f2,0.35,[.5 .6]);   colormap(gray);    
    set(gcf,'Color','white'); close(f1); close(f2);
    set(gcf, 'Position', [0 0 800 800]);	
    if(saveFigs)
        SaveFigure(['Density_eta_',strEta],v2struct(optsPhys,optsNum));
    end

end