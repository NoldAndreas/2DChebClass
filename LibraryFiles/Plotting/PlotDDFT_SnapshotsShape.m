function PlotDDFT_SnapshotsShape(input,file_name,opts)

    if(nargin < 3)
        opts = {}; %possible opts: 4Snapshots        
    end

    if(~isstr(input) && isa(input,'ContactLineHS'))        
        subArea    = input.subArea;
        optsPhys   = input.optsPhys;
        optsNum    = input.optsNum;
        data       = input.dynamicsResult;
        filename   = input.FilenameDyn;
        shape      = input.IDC;
    else
        v2struct(input);
    end
    
    %*****************************
    %Initialization of Data         
    v2struct(data);
    
    % overdamped case
    if(~exist('UV_t','var'))
        UV_t = flux_t;
    end
    
    plotTimes  = t;
    n          = length(plotTimes);
    fl_norm    = 0.1*max(max(max(abs(UV_t))));
    disp(['Flux normalized with ',num2str(fl_norm)]);
    
    %Choose times for snapshots
    if(IsOption(opts,'4Snapshots'))
        noPlots = 4;                
        cols    = 2;        
    else
        noPlots = 12;                                        
        cols    = 4;        
    end 
    ct      = 1:n;
    d       = floor((n-1)/(noPlots-1));                
	mark    = (mod(ct-1,d) == 0);
	rows    = ceil(sum(mark)/cols);        
    
	if(size(rho_t,3) == 1)
        rho_s(:,1,:)   = rho_t;  rho_t = rho_s;
        flux_s(:,1,:)  = UV_t;   flux_t = flux_s;        
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
    if(~IsOption(opts,'noNewFigure'))       
        figure('Color','white','Position', [0 0 (200+300*cols) (200+300*rows)]);            
    end
    %*****************************    
    
    pl_j = 1;
	for j=1:d:length(plotTimes)       
        
        rho     = rho_t(:,:,j);
        t       = plotTimes(j);
        
        subplot(rows,cols,pl_j);
        
        rhoMin = min(min(min(rho)));
        rhoMax = max(max(max(rho)));
        
        for iSpecies=1:nSpecies
             optsPlot.linecolor = lineColour{iSpecies};              
             
             shape.plotFlux(UV_t(:,iSpecies,j),shape.Ind.bound,fl_norm,1.5,lineColour{iSpecies}); hold on;                  
             shape.plotFlux(UV_t(:,iSpecies,j),~shape.Ind.bound,fl_norm,0.5,lineColour{iSpecies}); hold on;           
             
             if(isfield(optsPhys,'rhoLiq_sat') && isfield(optsPhys,'rhoGas_sat'))
                %optDetails.y2CartShift = -0.5;
                optDetails.clabel      = false;  
                optDetails.linewidth   = 1.5;  

                %optDetails.nContours = [0.1,0.2,0.3,0.4,0.5,0.6,0.7];        
                rhoLiq_sat = optsPhys.rhoLiq_sat;
                rhoGas_sat = optsPhys.rhoGas_sat;
                drho = rhoLiq_sat - rhoGas_sat;

                optDetails.nContours = rhoGas_sat + 0.1*drho;
                optDetails.linecolor = 'b';
                optDetails.linestyle = '--';
                shape.plot(rho,'contour',optDetails);  hold on;  

                optDetails.nContours = rhoGas_sat + 0.5*drho;
                optDetails.linecolor = [0 0.75 0];
                shape.plot(rho,'contour',optDetails);  hold on;  

                optDetails.nContours = rhoGas_sat + 0.9*drho;
                optDetails.linecolor = 'r';
                shape.plot(rho,'contour',optDetails);  hold on;                  
             else
                 optsPlot.nContours = [0.1,0.3,0.6];%rhoMin + [0.25,0.5,0.75]*(rhoMax - rhoMin);
                 shape.plot(rho(:,iSpecies),'contour',optsPlot); hold on;            
             end                          
             
             if(bool_subSp)
                subArea.PlotBorders();
             end

        end
        if(IsOption(opts,'NumericsManuscript'))                        
            xlabel('$y_1$','Interpreter','Latex','fontsize',20);
            ylabel('$y_2$','Interpreter','Latex','fontsize',20);
            title(['t = ', num2str(t,2)]);%,...%,'\\ j_{max} = ',num2str(max(max(abs(flux_t(:,:,j)))),'%.1e'),'$'],...            
        else
            title(['t = ', num2str(t,3)]);%,...%,'\\ j_{max} = ',num2str(max(max(abs(flux_t(:,:,j)))),'%.1e'),'$'],...        
        end
        %    'fontsize',15);
        pl_j = pl_j +1;                
    end        
    if((nargin > 1) && ~isempty(file_name))
        SaveFigure(file_name,v2struct(optsPhys,optsNum));
	end
end