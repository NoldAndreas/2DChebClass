function PlotDensitySlicesMovie(this)   

    %% Initialization 
    n     = 100;
    y2Max = 15;
    y1Min = 0;
    y1Max = 15;
    
    rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
    rhoGas_sat     = this.optsPhys.rhoGas_sat;     
    
    y1P = y1Min + (y1Max-y1Min)*(0:1:n)/n;
    this.optsNum.PlotArea = struct('y1Min',y1Min,'y1Max',y1Max,...
                                   'y2Min',0.5,'y2Max',y2Max,...
                                   'N1',100,'N2',100);
    InitInterpolation(this,true);


    k = 1; fileNames = [];

    %% Plotting
    figure('color','white','Position',[0 0 1200 600]);    

    for i = 1:n
        % get adsorption

        ell      = this.HS.doIntFLine([y1P(i) y1P(i)],[0.5 y2Max],this.rho_eq-rhoGas_sat,'CHEB')/(rhoLiq_sat-rhoGas_sat);
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);        

        subplot(1,2,1);
        hold off;
        plot(this.AdsorptionIsotherm_Pts.y2_kv-0.5,rho,'b','linewidth',1.5); hold on;
        this.HS.doPlotFLine([y1P(i) y1P(i)],[0.5 y2Max],this.rho_eq,'CART',true);
        xlim([0 (y2Max-0.5)]);
        ylim([0 1.1]);%max(this.rho_eq)]);
        xlabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
        ylabel('$n\sigma^3$','Interpreter','Latex','fontsize',25);
        set(gca,'fontsize',20);                        
        set(gca,'linewidth',1.5);          
        pbaspect([(y1Max-y1Min) (y2Max) 1]);

        subplot(1,2,2);
        hold off;
        PlotEquilibriumResults(this,[y1Min y1Max],[0.5 y2Max],true,false); hold on;        
        plot([y1P(i) y1P(i)],[0 (y2Max-0.5)],'k--','linewidth',1.5);

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

    end

    %% Save Movie
    str         = 'DensitySlices';
    allPdfFiles = [str,'.pdf'];
    swfFile     = [str,'.swf'];

    system(['C:\pdftk.exe ', fileNames ,' cat output ',allPdfFiles]);    
    system(['C:\pdf2swf.exe -s framerate=5 -o ',swfFile,' ', allPdfFiles]);
    system(['copy ',getPDFMovieFile('Movie1',1),' ',str,'POSTER.pdf']);
    system(['del ',fileNames]);       
    disp(['Swf Movie` saved in: ',swfFile]);

end