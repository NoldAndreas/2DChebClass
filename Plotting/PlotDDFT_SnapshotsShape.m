function PlotDDFT_SnapshotsShape(input,file_name)

    %*****************************
    %Initialization of Data     
    v2struct(input);
    v2struct(data);
    
    plotTimes  = t;%optsNum.plotTimes;    
    fl_norm    = 0.1*max(max(max(abs(flux_t))));
    disp(['Flux normalized with ',num2str(fl_norm)]);
    
    %Choose times for snapshots
    noPlots = 12;
    n       = length(plotTimes);
    d       = ceil((n-1)/(noPlots-1));
    ct      = 1:n;
    mark    = (mod(ct-1,d) == 0);
    
    cols = 4;
    rows = ceil(sum(mark)/cols);
        
	if(size(rho_t,3) == 1)
        rho_s(:,1,:)   = rho_t;  rho_t = rho_s;
        flux_s(:,1,:)  = flux_t; flux_t = flux_s;        
    end
	nSpecies = size(rho_t,2);
    
	if(exist('optsPlot','var') && isstruct(optsPlot) && isfield(optsPlot,'lineColourDDFT'))
        lineColour=optsPlot.lineColourDDFT;
    else
        lineColour={'b','r','g','k','m'};
    end   
    
       if(isfield(data,'Subspace'))
        v2struct(Subspace);        
        IP            = shape.SubShapePts(subArea.Pts);
        Int_SubOnFull = subArea.ComputeIntegrationVector()*IP;
        bool_subSp = true;
    else
        bool_subSp = false;
    end
    %**************************************
    %Initialization of figure and screen, and movie
    close all;
    figure('Color','white','Position', [0 0 1500 (200+300*rows)]);            
    %*****************************    
    
    pl_j = 1;
	for j=1:d:length(plotTimes)       
        
        rho     = rho_t(:,:,j);
        t       = plotTimes(j);
        
        subplot(rows,cols,pl_j);
        
        for iSpecies=1:nSpecies
             optsPlot.linecolor = lineColour{iSpecies};              
             
             shape.plotFlux(flux_t(:,iSpecies,j),shape.Ind.bound,fl_norm,1.5,lineColour{iSpecies}); hold on;                  
             shape.plotFlux(flux_t(:,iSpecies,j),~shape.Ind.bound,fl_norm,0.5,lineColour{iSpecies}); hold on;           
             shape.plot(rho(:,iSpecies),'contour',optsPlot); hold on;            
             if(bool_subSp)
                subArea.PlotBorders();
             end

        end
        title(['t = ', num2str(t)]);%,...%,'\\ j_{max} = ',num2str(max(max(abs(flux_t(:,:,j)))),'%.1e'),'$'],...
        %    'fontsize',15);
        pl_j = pl_j +1;                
    end        
    if(nargin > 1)
        SaveFigure(file_name);
            %print2eps(file_name,gcf);
            %saveas(gcf,[file_name,'.fig']);
            %disp(['Snapshots saved in ',file_name,'.eps/.fig.']);
	end
end