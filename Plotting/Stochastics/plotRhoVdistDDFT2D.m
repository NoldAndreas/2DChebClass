function handles=plotRhoVdistDDFT2D(rho,v,IDC,optsPlot,plotHandles)
%plotRhoVdistDDFT(rhov,r,optsPlot,plotHandles,type)
%   plots position and velocity distributions in axes given by plotHandles
%
% INPUTS: 
%  rhov         -- (2*N,1) vector with columns [rho; v] specified  
%                   at grid points r
%  r            -- (N,1) vector of grid points
%  optsPlot     -- a structure of size [1 1] containing:
%                        lineStyle  (linestyles for stochastic plots, 
%                                    should be  of the form ['b  ';'--b'])
%                        plotDensity (true/false whether to plot the
%                                        density; false->distruibution)
%                        geom
%                        N1
%                        N2
%  plotHandles  -- a structure of size [1 1] containing:
%                        hRa        (handle for position axis)
%                        hPa        (handle for momentum axis)
%  type         -- either 'r' or 'rv' determining if we do momentum plots
%
% OUTPUTS:
%  handles      -- a structure containing the handles of the plot lines,
%                  hr for the position and if type=rv also hp for velocity

%plotRhoVdistDDFT2D(rhoEnd,[],temp.Interp.pts1,temp.Interp.pts1,opts,[],'r')

Interp = IDC.Interp;

PlotArea.y1Min = max(optsPlot.rMin(1),min(IDC.Pts.y1));
PlotArea.y2Min = max(optsPlot.rMin(2),min(IDC.Pts.y2));
PlotArea.y1Max = min(optsPlot.rMax(1),max(IDC.Pts.y1));
PlotArea.y2Max = min(optsPlot.rMax(2),max(IDC.Pts.y2));
PlotArea.N1 = 10;
PlotArea.N2 = 10;

InterpFlux = IDC.InterpolationPlotCart(PlotArea,false);

contourWidth=get(0,'defaultlinelinewidth');

nSpecies=size(rho,2);

Nplot1=Interp.Nplot1;
Nplot2=Interp.Nplot2;

%type=optsPlot.type;

% if no figure/axes, create them
if(isempty(plotHandles))
     
    fullscreen = get(0,'ScreenSize');

    hRPf=figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
    % set background colour to white
    set(hRPf,'Color','w');

    switch optsPlot.plotType
        case 'surf'
            handles=tightsubplot(2,1,0.075,0.075,0.075);
            hRa=handles(1)*ones(1,nSpecies);
            hPa=handles(2)*ones(1,nSpecies);
        case 'contour'
            hRa=axes;
            % invisible axes off the figure
            hPa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
    end

    % figure and axes handles to be passed to plotting functions
    handles=struct('hRPf',hRPf,'hRa',hRa,'hPa',hPa);
      
    % otherwise assign axis handles
else
    hRa=plotHandles.hRa;
    hPa=plotHandles.hPa;
    
    if(length(hRa)<nSpecies)
        for iSpecies=1:nSpecies
            hRa(iSpecies) = hRa(1);
        end
    end

    if(length(hPa)<nSpecies)
        for iSpecies=1:nSpecies
            hPa(iSpecies) = hPa(1);
        end
    end
    
end


faceColour=optsPlot.faceColour;
lineStyle = optsPlot.lineStyle;

x1=Interp.pts1;
x2=Interp.pts2;

x1f = InterpFlux.pts1;
x2f = InterpFlux.pts2;


switch optsPlot.geom

    case 'polar2D'
        % x1=r, x2=theta
        y1 = x1.*cos(x2); 
        y2 = x1.*sin(x2);
        
        y1f = x1f.*cos(x2f); 
        y2f = x1f.*sin(x2f);
    case 'planar2D'
        y1=x1;
        y2=x2;
        
        y1f=x1f;
        y2f=x2f;
        
end

y1=reshape(y1,Nplot2,Nplot1);
y2=reshape(y2,Nplot2,Nplot1);


N1=Interp.N1;
N2=Interp.N2;

hr=zeros(nSpecies,1);
hp=hr;

y1pad=optsPlot.rMin(1)-1;
y2pad=optsPlot.rMin(2)-1;

fluxNorm=optsPlot.fluxNorm;


for iSpecies=1:nSpecies
    
    rhoS=rho(:,iSpecies);
    
    
    fluxS=v(:,iSpecies);
    fluxS1temp = fluxS(1:N1*N2);
    fluxS2temp = fluxS(N1*N2+1:end);
    
    switch optsPlot.geom

        case 'polar2D'
            % flux is given in cartesian even in polar case

            fluxS1=fluxS1temp;
            fluxS2=fluxS2temp;

        case 'planar2D'
        
            fluxS1=fluxS1temp;
            fluxS2=fluxS2temp;
    end
 
    fluxS1interp = InterpFlux.InterPol*fluxS1;
    fluxS2interp = InterpFlux.InterPol*fluxS2;
 
    rhoS = real(Interp.InterPol*rhoS);
    
    rhoS=reshape(rhoS,Nplot2,Nplot1);
    
    
    switch optsPlot.plotType
        case 'surf'
            % surface plot

            hr(iSpecies)=surf(hRa(iSpecies),y1,y2,rhoS); 
            set(hr(iSpecies),'FaceColor',faceColour{iSpecies});
            alpha(hr(iSpecies),0.5);
            
            hCa=hPa(iSpecies);
        
        case 'contour'
            
            hCa=hRa(iSpecies);
    end
    
    cutoff=10^(-10);
    
    rhoS(rhoS<cutoff)=0;
    
    [C,h]=contour(hCa,y1,y2,rhoS);     
    
    if(strcmp(lineStyle{iSpecies},'none'))
        lineStyle{iSpecies} = '-';
    end
    
    
    set(h,'color',faceColour{iSpecies},'linewidth',contourWidth,'linestyle',lineStyle{iSpecies});
    clabel(C,h,'Color',faceColour{iSpecies});
    set(h,'ShowText','off');
    hold(hCa,'on');
    
    %quiver(ha,[y1f(mask)],[y2f(mask)],[fluxS1(mask)],[fluxS2(mask)]); 
    
    %h=quiver(hCa,[y1pad ; y1f(mask)],[y2pad; y2f(mask)],[fluxNorm; fluxS1(mask)],[0; fluxS2(mask)]); 
    
    %h=quiver(hCa,[y1pad ; y1f(mask)],[y2pad; y2f(mask)],[fluxNorm; fluxS1interp(mask)],[0; fluxS2interp(mask)]); 
    
    if(isfield(optsPlot,'doFlux') && ~optsPlot.doFlux)
        % don't plot flux
    else
        h=quiver(hCa,[y1pad ; y1f],[y2pad; y2f],[fluxNorm; fluxS1interp],[0; fluxS2interp]); 
    end
    
    set(h,'color',faceColour{iSpecies});
    axis(hCa,'equal');    
       
    hold(hRa(iSpecies),'on')
    hold(hPa(iSpecies),'on')
    
end

% set output handles of surfaces
handles.hr=hr;
handles.hp=hp;

