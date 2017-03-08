function PlotDDFT1D(input,Bool_Record)

    v2struct(input);
    v2struct(data);

    plotTimes = optsNum.plotTimes;

    if(size(rho_t,3)==1)
        rho_s(:,1,:)   = rho_t;  rho_t = rho_s;
        flux_s(:,1,:)  = flux_t; flux_t = flux_s;
    end

    nSpecies = size(rho_t,2);

    if(nSpecies>1)
        nRows    = ceil(nSpecies)+1;
    else
        nRows    = 3;
    end

    rho_ic = rho_t(:,:,1);

    if(exist('optsPlot','var') && isstruct(optsPlot) && isfield(optsPlot,'lineColourDDFT'))
        lineColour=optsPlot.lineColourDDFT;
    else
        lineColour={'r','b','g','k','m'};
    end

    if(nargin == 1)
        Bool_Record = false;
    end

    figure

    for i=1:length(plotTimes)
        rho     = rho_t(:,:,i);
        flux    = flux_t(:,:,i);
        t       = plotTimes(i);

        h1 = subplot(nRows,2,1);
        cla(h1);
        optsPlot.linestyle = '-';
        optsPlot.dist = true;
        for iSpecies=1:nSpecies
             optsPlot.linecolor = lineColour{iSpecies};
             shape.plot(rho(:,iSpecies),optsPlot); hold on;
        end


%       Check mass conservation
        subplot(nRows,2,2)
         for iSpecies=1:nSpecies
            rho_diff = (rho(:,iSpecies)-rho_ic(:,iSpecies));
            plot(t,shape.Int*rho_diff,'o','Color',lineColour{iSpecies});
            hold on;
         end

         xlabel('t');
         title(['Mass in the full system: ',num2str(shape.Int*rho)]);


         for iSpecies=1:nSpecies
            optsPlot.linestyle = '-';
            optsPlot.linecolor = lineColour{iSpecies};
            optsPlot.dist = true;
            subplot(nRows,2,2+iSpecies+2*(ceil(iSpecies/2)-1));
            shape.plot(rho(:,iSpecies),optsPlot);
            hold off
            title(['Species ' num2str(iSpecies)]);


            subplot(nRows,2,2+iSpecies+2*ceil(iSpecies/2));
            optsPlot.linestyle = '-';
            optsPlot.dist = false;
            shape.plot(flux(:,iSpecies),optsPlot);
            hold off
            title(['Flux ' num2str(iSpecies)]);
         end


         pause(0.05);

    end

end