function PlotDDFT(input,Bool_Record)

    %*****************************
    %Initialization of Data     
    v2struct(input);
    v2struct(data);
    
%     if(isfield(optsNum.PhysArea,'N'))
%         N1 = optsNum.PhysArea.N(1);
%         N2 = optsNum.PhysArea.N(2);
%     else
%         N1 = optsNum.PhysArea.N1; N2 = optsNum.PhysArea.N2;
%     end
    plotTimes  = optsNum.plotTimes;          
    
    if(isfield(data,'Subspace'))
        v2struct(Subspace);        
        IP            = shape.SubShapePts(subArea.Pts);
        Int_SubOnFull = subArea.ComputeIntegrationVector()*IP;
        bool_subSp = true;
    else
        bool_subSp = false;
    end
    
    if(size(rho_t,3) == 1)
        %X_s(:,1,:)   = X_t;  X_t = X_s;
        rho_s(:,1,:)   = rho_t;  rho_t = rho_s;
        flux_s(:,1,:)  = flux_t; flux_t = flux_s;        
    end
    nSpecies = size(rho_t,2);
    nRows    = ceil(nSpecies/2)+1;
           
    rho_ic = rho_t(:,:,1);

    fl_norm = 0.1*max(max(max(abs(flux_t))));    
    
    if(exist('optsPlot','var') && isstruct(optsPlot) && isfield(optsPlot,'lineColourDDFT'))
        lineColour=optsPlot.lineColourDDFT;
    else
        lineColour={'r','b','g','k','m'};
    end
    cornerColour=cell(nSpecies,1);
    for iSpecies=1:nSpecies
        cornerColour{iSpecies}='m';
    end        
    
    %**************************************
    %Initialization of figure and screen, and movie
    if(nargin == 1)
        Bool_Record = false;
    end

    close all;
    figure       
            
    %*****************************
    if((islogical(Bool_Record) && Bool_Record) || ischar(Bool_Record))
        if(ischar(Bool_Record))
            gifFile = [Bool_Record,'.gif'];
        else
            gifFile = getMovieFile('Movie');    
        end
        
        %(a) Plot Snapshots
        PlotDDFT_SnapshotsShape(input,[gifFile(1:end-4) '_Snapshots']);        
        
        %(b1) Plot Mass evolution
        figure('Color','white','Position', [0 0 800 800]); %Set background color        
        rho     = permute(rho_t(:,1,:),[1 3 2]);
        for iSpecies=1:nSpecies            
            rho_diff = rho-rho_ic(:,iSpecies)*ones(1,length(plotTimes));
            plot(plotTimes,shape.Int*rho_diff,'Color',lineColour{iSpecies},'linewidth',1.5); hold on; 
            %plot(plotTimes,shape.Int*rho_diff,'o','Color',lineColour{iSpecies}); hold on; 
            if(bool_subSp)
                plot(plotTimes,Int_SubOnFull*rho_diff+accFlux','m','linewidth',1.5);
                legend('Full Domain','Subdomain','location','NorthWest');
            end
        end
         grid on;
         xlabel('t','fontsize',20);
         ylabel('Mass','fontsize',20);         
         set(gca,'fontsize',20);                        
         set(gca,'linewidth',1.5);      
         
         print2eps([Bool_Record , '_Mass'],gcf);
         saveas(gcf,[Bool_Record , '_Mass.fig']);
         
        %(b2) Plot Mass in Subbox
        if(bool_subSp)
            figure('Color','white','Position', [0 0 800 800]); %Set background color        
            rho     = permute(rho_t(:,1,:),[1 3 2]);
            for iSpecies=1:nSpecies                 
                plot(plotTimes,Int_SubOnFull*rho,'m','linewidth',1.5);                
            end
            grid on;
            xlabel('t','fontsize',20);
            ylabel('Mass','fontsize',20);         
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);      

            print2eps([Bool_Record , '_Mass_SubArea'],gcf);
            saveas(gcf,[Bool_Record , '_Mass_SubArea.fig']); 
        end
         
        %(c) Plot Movie        
        figure('Color','white','Position', [0 0 1500 1000]); %Set background color        
        k = 1; fileNames = [];
        dT = ceil(length(plotTimes)/20);
        T_n_Max = length(plotTimes);
        for i=1:dT:T_n_Max
            
            t       = plotTimes(i);

            %Plot Interpolation
            hold off;
            %shape.doPlots(rho,'contour');  hold on;
            %shape.doPlotsFlux(flux_t(:,1,i),~shape.Ind.bound,fl_norm,0.5,'k'); hold on;            
            for iSpecies=1:nSpecies    
                rho     = rho_t(:,iSpecies,i);
                optsPlot.linecolor = lineColour{iSpecies}; 
                shape.doPlotsFlux(flux_t(:,iSpecies,i),~shape.Ind.bound,fl_norm,0.5,lineColour{iSpecies}); hold on;           
                shape.doPlots(rho,'contour',optsPlot); hold on;                         
            end
            %plot([0;0],[i;i]);
            title(['t = ', num2str(round(t))]);               
            set(gca,'fontsize',30);
            set(gca,'linewidth',3);
            %view([2,5,2]);
            
            h = get(gca,'xlabel'); set(h,'fontsize',35);
            h = get(gca,'ylabel'); set(h,'fontsize',35);
            h = get(gca,'title');  set(h,'fontsize',35);
        
            %For gif-recording
            Record(i,gifFile);
            
            %for swf Recording
            fileName = getPDFMovieFile('Movie1',k);
            %fileName = ['Movie1',num2str(k),'.pdf'];
            save2pdf(fileName,gcf);
            k = k+1;
            fileNames = [fileNames,' ',fileName];
        end
        
        disp(['Gif Movie` saved in: ',gifFile]);
        
        str         = Bool_Record;
        allPdfFiles = [str,'.pdf'];
        swfFile     = [str,'.swf'];

        system(['C:\pdftk.exe ', fileNames ,' cat output ',allPdfFiles]);    
        system(['C:\pdf2swf.exe -s framerate=5 -o ',swfFile,' ', allPdfFiles]);
        system(['copy ',getPDFMovieFile('Movie1',1),' ',str,'POSTER.pdf']);
        system(['del ',fileNames]);       
        disp(['Swf Movie` saved in: ',swfFile]);
        
        return;
    end           
    
    dT = ceil(length(plotTimes)/100);
    for i=1:dT:length(plotTimes)                
        rho     = rho_t(:,:,i);
        t       = plotTimes(i);
                
        %Plot Fluxes
        h1 = subplot(nRows,2,1);    
        cla(h1);
        for iSpecies=1:nSpecies
             optsPlot.linecolor = lineColour{iSpecies}; 
             shape.doPlotsFlux(flux_t(:,iSpecies,i),shape.Ind.bound,fl_norm,1.5,lineColour{iSpecies}); hold on;                  
             shape.doPlotsFlux(flux_t(:,iSpecies,i),~shape.Ind.bound,fl_norm,0.5,lineColour{iSpecies}); hold on;           
             shape.doPlots(rho(:,iSpecies),'contour',optsPlot); hold on;            
             if(bool_subSp)
                subArea.PlotBorders();
             end                
        end
        title(['t = ', num2str(t),' max(flux) = ',num2str(max(max(abs(flux_t(:,:,i)))))]);
        
%       Check mass conservation
        subplot(nRows,2,2)   
         hold on;
         for iSpecies=1:nSpecies
            rho_diff = (rho(:,iSpecies)-rho_ic(:,iSpecies));
            plot(t,shape.Int*rho_diff,'o','Color',lineColour{iSpecies}); hold on;                                    
            if(bool_subSp)
                plot(t,Int_SubOnFull*rho_diff+accFlux(i),'om');
                legend('Full Domain','Subdomain');
            end
         end
         xlabel('t');
         title(['Mass in the full system: ',num2str(shape.Int*rho)]);                
         
 %         Plot Density profiles         
         hold off;
         for iSpecies=1:nSpecies
            subplot(nRows,2,2+iSpecies);
            shape.doPlots(rho(:,iSpecies),'SC');
            view([-2,5,2]);     
            title(['Species ' num2str(iSpecies)]);
         end

         pause(0.05);        

        
%         
%         maxBoundFlux = max(fl_y1(maB).^2 + fl_y2(maB).^2);
%         maxCornerFlux = max(fl_y1(maC).^2 + fl_y2(maC).^2);
%         if( maxCornerFlux > 10^(-10)) 
%             NormQuiverPlot(Pts,fl,maC,[],[],2,'m');
%             %quiver(y1_s(maC),y2_s(maC),fl_y1(maC),fl_y2(maC),'LineWidth',2,'Color','m'); hold on;        
%         end
%         
%         doPlots_IP_Contour(data.Interp,rho);        hold on;         
%          xlim(xl); ylim(yl);
%          xlabel('y1');        ylabel('y2');
%          
%         title(['Flux, max [all,boundary,corner]: ' , num2str(max(fl_y1.^2 + fl_y2.^2)) ,...
%                ' , ' , num2str(maxBoundFlux),...
%                ' , ' , num2str(maxCornerFlux)]);
% 
%         
%         if(isfield(data,'Subspace'))
%             v2struct(Subspace);
%             
%         Pathy1 = [Path(1).pts_y1;Path(2).pts_y1;Path(3).pts_y1;Path(4).pts_y1];
%         Pathy2 = [Path(1).pts_y2;Path(2).pts_y2;Path(3).pts_y2;Path(4).pts_y2];
% 
%         IPFlux1 =  [Path(1).InterpolFlux(1:end/2,:);Path(2).InterpolFlux(1:end/2,:);...
%                     Path(3).InterpolFlux(1:end/2,:);Path(4).InterpolFlux(1:end/2,:)];
%         IPFlux2 =  [Path(1).InterpolFlux(end/2+1:end,:);Path(2).InterpolFlux(end/2+1:end,:);...
%                     Path(3).InterpolFlux(end/2+1:end,:);Path(4).InterpolFlux(end/2+1:end,:)];
% 
%         maS = (Pathy1 < max(y1Plot));
%         quiver(Pathy1(maS),Pathy2(maS),IPFlux1(maS,:)*fl,IPFlux2(maS,:)*fl,'LineWidth',2,'Color','green');  hold off;
%         
%         %Plot Flux through boundary
%         subplot(2,2,4);   
%         hold off;
%         %Right boundary
%         plot(InterpPath(1).pts_x,InterpPath(1).InterPol*fl,'');  hold on;
%         plot(P``ath(1).pts_x,Path(1).InterPol*fl,'o'); hold on;
%                 
%         %Bottom boundary
%         plot(InterpPath(2).pts_x,InterpPath(2).InterPol*fl,':');  hold on;
%         plot(Path(2).pts_x,Path(2).InterPol*fl,'o');
%         
%         %Upper boundary
%         plot(InterpPath(3).pts_x,InterpPath(3).InterPol*fl,'-.');  hold on;
%         plot(Path(3).pts_x,Path(3).InterPol*fl,'o');
%         
%         %Left boundary
%         plot(InterpPath(4).pts_x,InterpPath(4).InterPol*fl,'--');  hold on;
%         plot(Path(4).pts_x,Path(4).InterPol*fl,'o');
%         
%         legend('Right','','Bottom','','Upper','','Left','');
%         title(['Flux through the path is: ' , num2str(Int_of_path*fl)]);
%         xlabel('y');
%         end
        
        pause(0.02);        
    end          
    
end