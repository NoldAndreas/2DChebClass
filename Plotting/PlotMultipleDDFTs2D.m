function PlotMultipleDDFTs2D(struct)

    nDDFT = length(struct);

    plotTimes = struct(1).optsNum.plotTimes;

    if(nDDFT>1)
     for iDDFT = 2:nDDFT
         plotTimesI = struct(iDDFT).optsNum.plotTimes;
         if(~isequal(plotTimesI,plotTimes))
             fprintf(1,'plotTimes not equal\n');
             return;
         end
     end
    end
   
    nSpeciesMax = 1;
     
    for iDDFT = 1:nDDFT
        rho_temp = struct(iDDFT).data.rho_t;


        nSpeciesMax = max(nSpeciesMax,size(rho_temp,2));

        struct(iDDFT).nSpecies = size(rho_temp,2);

        if(~isfield(struct(iDDFT),'optsPlot') ...
                || (isfield(struct(iDDFT),'optsPlot') ...
                        && ~isfield(struct(iDDFT).optsPlot,'lineColourDDFT') ) )
            struct(iDDFT).optsPlot.lineColourDDFT={'r','b','g','k','m'};
        end

    end
         
    
%     if(nSpeciesMax>1)
%         nRows    = ceil(nSpeciesMax)+1;
%     else
%         nRows    = 3;
%     end
      

    nRows = 1;

    figure
    h1 = subplot(nRows,2,1); 
    h2 = subplot(nRows,2,2); 
    
    fl_norm = 0;
    for iDDFT = 1:nDDFT
        flux = struct(iDDFT).data.flux_t;
        
        fl_norm = max(0,0.1*max(max(max(abs(flux)))));
    end

    for i=1:length(plotTimes)
        
        cla(h1);
        cla(h2);
%        cla(hV);
       
        for iDDFT = 1:nDDFT
            
            rho     = struct(iDDFT).data.rho_t(:,:,i);
            rho_ic  = struct(iDDFT).data.rho_t(:,:,1);
            flux    = struct(iDDFT).data.flux_t(:,:,i);
            %V       = struct(iDDFT).data.V_t(:,:,i);

            t       = plotTimes(i);
            
            axes(h1);  
            struct(iDDFT).optsPlot.linestyle = '-';

            %Plot densities
            for iSpecies=1:struct(iDDFT).nSpecies
                struct(iDDFT).optsPlot.linecolor = struct(iDDFT).optsPlot.lineColourDDFT{iSpecies};
                struct(iDDFT).data.shape.plot(rho(:,iSpecies),'SC',struct(iDDFT).optsPlot); 
                hold on;     
            end
            if(iDDFT==nDDFT)
                hold off;
            end
            
            %Plot fluxes
            axes(h2);
            for iSpecies=1:struct(iDDFT).nSpecies
                struct(iDDFT).optsPlot.linecolor = struct(iDDFT).optsPlot.lineColourDDFT{iSpecies};
                struct(iDDFT).data.shape.plotFlux(flux(:,iSpecies),struct(iDDFT).data.shape.Ind.bound,fl_norm,1.5, ...
                            struct(iDDFT).optsPlot.lineColourDDFT{iSpecies});
                hold on;               
                struct(iDDFT).data.shape.plotFlux(flux(:,iSpecies),~struct(iDDFT).data.shape.Ind.bound,fl_norm,0.5, ...
                            struct(iDDFT).optsPlot.lineColourDDFT{iSpecies});
                hold on;
                struct(iDDFT).data.shape.plot(rho(:,iSpecies),'contour', struct(iDDFT).optsPlot);
                hold on;     
            end
            if(iDDFT==nDDFT)
                hold off;
            end
            
             

            %Plot Density profiles         
%             for iSpecies=1:struct(iDDFT).nSpecies
%                 struct(iDDFT).optsPlot.linestyle = '-';
%                 struct(iDDFT).optsPlot.linecolor = struct(iDDFT).optsPlot.lineColourDDFT{iSpecies};
%                 subplot(nRows,2,2+iSpecies+2*(ceil(iSpecies/2)-1));
%                 struct(iDDFT).data.shape.plot(rho(:,iSpecies),struct(iDDFT).optsPlot);
%                 hold on;
%                 if(iDDFT==nDDFT)
%                     hold off;
%                 end
%                 title(['Species ' num2str(iSpecies)]);
%             end
%             
% 
%             
%             for iSpecies=1:struct(iDDFT).nSpecies
%                 subplot(nRows,2,2+iSpecies+2*ceil(iSpecies/2));
%                 struct(iDDFT).optsPlot.linestyle = '-';
%                 struct(iDDFT).optsPlot.linecolor = struct(iDDFT).optsPlot.lineColourDDFT{iSpecies};
%                 struct(iDDFT).optsPlot.dist = false;
%                 struct(iDDFT).data.shape.plot(flux(:,iSpecies),struct(iDDFT).optsPlot);
%                 hold on
%                 if(iDDFT==nDDFT)
%                     hold off;
%                 end
%                 title(['Flux ' num2str(iSpecies)]);
%             end            
           
        end
        %pause         
        pause(0.05);                 
    end 
end

