function PlotDensitySlicesMovie(this)   

    global dirData
    %% Initialization 
    n     = 50;
    y2Max = 15;
    y1Min = min(this.y1_SpectralLine.Pts.y);
    y1Max = max(this.y1_SpectralLine.Pts.y);
    
    rho_eq = this.GetRhoEq;
    
    rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
    rhoGas_sat     = this.optsPhys.rhoGas_sat;     
    mu_sat         = this.optsPhys.mu_sat;     
    
    y1P = y1Min + (y1Max-y1Min)*(0:1:n)/n;
    this.optsNum.PlotArea = struct('y1Min',y1Min,'y1Max',y1Max,...
                                   'y2Min',0.5,'y2Max',y2Max,...
                                   'N1',100,'N2',100);
    %InitInterpolation(this);
    
    if(IsDrying(this))
        rho_bulk = rhoLiq_sat;
    else
        rho_bulk = rhoGas_sat;
    end
    
    k = 1; fileNames = [];

    %% Plotting

    for i = 1:n
        % get adsorption

        ell      = this.IDC.doIntFLine([y1P(i) y1P(i)],[0.5 y2Max],rho_eq-rho_bulk,'CHEB')/(rhoLiq_sat-rhoGas_sat);
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);        
        
        f1 = figure('Position',[0 0 800 800]);
        plot(this.AdsorptionIsotherm.mu-mu_sat,this.AdsorptionIsotherm.FT,'b','linewidth',1.5); hold on;
        plot(mu-mu_sat,ell,'or','MarkerSize',7,'MarkerFace','r');
        set(gca,'fontsize',15);                        
        set(gca,'linewidth',1.5);          
        xlabel('$\Delta \mu/\varepsilon$','Interpreter','Latex','fontsize',18);
        ylabel('$\ell/\sigma$','Interpreter','Latex','fontsize',18);
        
        xlim([(min(this.AdsorptionIsotherm.mu)-0.02) (max(this.AdsorptionIsotherm.mu)+0.02)]-mu_sat);        
        xTic = [(min(this.AdsorptionIsotherm.mu)) (max(this.AdsorptionIsotherm.mu))]-mu_sat;
        set(gca,'XTick',round(xTic*100)/100);        
        
        f2 = figure('color','white','Position',[0 0 1200 600]);   
        
        subplot(1,2,1);
        hold off;
        plot(this.AdsorptionIsotherm.Pts.y2-0.5,rho,'b','linewidth',1.5); hold on;        
        this.IDC.plotLine([y1P(i) y1P(i)],[0.5 y2Max],rho_eq,struct('dist0',true,'plain',true));
        xlim([0 (y2Max-0.5)]);
        ylim([0 1.1]);%max(this.rho_eq)]);
        xlabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
        ylabel('$n\sigma^3$','Interpreter','Latex','fontsize',25);
        set(gca,'fontsize',20);                        
        set(gca,'linewidth',1.5);          
        pbaspect([(y1Max-y1Min) (y2Max) 1]);
        pbaspect([1 1 1]);
        
        
        %f3 =  figure('color','white','Position',[0 0 1200 600]);    
        %s1 = subplot(1,2,1);
        
        %fig1 = get(ax2,'children');
        %copyobj(fig1,s1);

        subplot(1,2,2);
        hold off;
        PlotContourResults(this,true); hold on;        
        %PlotEquilibriumResults(this,[y1Min y1Max],[0.5 y2Max],true,false); hold on;        
        plot([y1P(i) y1P(i)],[0 (y2Max-0.5)],'k--','linewidth',1.5);
        
        
        inset2(f2,f1,[0.17,0.32],[0.34,0.6]);
        close(f1);

%         subplot(2,2,3);
%         hold off;
%         plot(this.AdsorptionIsotherm_Mu,this.AdsorptionIsotherm_FT,'b','linewidth',1.5); hold on;
%         plot(mu,ell,'o','MarkerEdgeColor','r','MarkerFaceColor','r');
%         
%         subplot(2,2,4);
%         hold off;
%         plot();

      %  pause(0.1);
        %For gif-recording
       % Record(i,gifFile);

        % for swf Recording
        fileName = getPDFMovieFile('Movie1',k);
        %fileName = ['Movie1',num2str(k),'.pdf'];
        save2pdf(fileName,gcf);
        k = k+1;
        fileNames = [fileNames,' ',fileName];           
        
        close all;
    end

    %% Save Movie
    str         = [dirData filesep 'Equilibrium' filesep this.FilenameEq ,'_DensitySlices'];
    allPdfFiles = [str,'.pdf'];
    swfFile     = [str,'.swf'];

    system(['C:\pdftk.exe ', fileNames ,' cat output ',allPdfFiles]);    
    system(['C:\pdf2swf.exe -s framerate=5 -o ',swfFile,' ', allPdfFiles]);
    system(['copy ',getPDFMovieFile('Movie1',1),' ',str,'POSTER.pdf']);
    system(['del ',fileNames]);       
    disp(['Swf Movie` saved in: ',swfFile]);

end