function outputFiles = plotMassChange2D(stoc,ddft,optsPlot)

fprintf(1,'Making mass change plots ... ');

% number of stochastic and DDFT plots to do
nStoc=size(stoc,2);
nDDFT=size(ddft,2);

% easy way to check if any of the computations include momentum
doP=~isempty(optsPlot.legTextP);

% fullscreen = get(0,'screensize');
% 
% hRPf=figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
% % set background colour to white
% set(hRPf,'Color','w');
% 
% hRa=axes;
% % invisible axes off the figure
% hPa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
% 
% % figure and axes handles to be passed to plotting functions
% handlesRP=struct('hRPf',hRPf,'hRa',hRa,'hPa',hPa); %#ok

if(nDDFT>0)
    lineColourDDFT=optsPlot.lineColourDDFT;
    lineStyleDDFT=optsPlot.lineStyleDDFT;
    % get type info -- used to decide whether to plot the velocity
    DDFTType=optsPlot.DDFTType;
else
    lineColourDDFT=[];
    lineStyleDDFT=[];
    % get type info -- used to decide whether to plot the velocity
    DDFTType=[];
end
    
% and file to save in
outputFile=[optsPlot.plotDir filesep 'massChange.pdf'];
    
%----------------------------------------------------------------------
% DDFT data plots
%----------------------------------------------------------------------

plotTimes = optsPlot.plotTimes;

subMasses = zeros(length(plotTimes),nDDFT);

hMf = figure('Position',[0,0,1200,800]);
hMa = axes;

for iDDFT=1:nDDFT

%     accFlux = DDFTStruct(1).dynamicsResult.Subspace.accFlux;
%     
%     plot(plotTimes,-cumsum(accFlux)*dt,'r')
    
    % get rho, v and r values
    rho=ddft(iDDFT).dynamicsResult.rho_t;

    optsPlot.type=DDFTType{iDDFT};

    IDC = ddft(iDDFT).IDC;
    subArea = ddft(iDDFT).subArea;
    InterPolSub = IDC.InterpolationMatrix_Pointwise(subArea.Pts.y1_kv,subArea.Pts.y2_kv);
    
    for iPlot = 1:length(plotTimes)

        
    
        % get values at appropriate time
        rhot=rho(:,:,iPlot);
        rhotSub = InterPolSub*rhot;

%         hold(hRa,'off');
        
%         if(strcmp(optsPlot.type,'rv'))
%             v=ddft(iDDFT).v_IP;
%             vt=v(:,:,plotPos(iPlot));
%         else
%             v=ddft(iDDFT).dynamicsResult.flux_t;
%             fluxNorm = 0.1*max(max(max(v)));
%             optsPlot.fluxNorm=fluxNorm;
% 
%             vt=v(:,:,1);
%         end
% 
%         optsPlot.faceColour=lineColourDDFT{iDDFT};
%         optsPlot.lineColour=lineColourDDFT{iDDFT};
%         optsPlot.lineStyle=lineStyleDDFT{iDDFT};

%         IDC = ddft(iDDFT).IDC;
% 
%         %plotRhoVdistDDFT2D(rhot,vt,IDC,optsPlot,handlesRP);
% 
%         hold(hRa,'on');
%         hold(hPa,'on');
% 
%         subArea = ddft(iDDFT).subArea;
%         InterPolSub = IDC.InterpolationMatrix_Pointwise(subArea.Pts.y1_kv,subArea.Pts.y2_kv);
% 
%         rhotSub = InterPolSub*rhot;
% 
%         v1 = vt(1:end/2);
%         v2 = vt(end/2+1:end);
%         v1Sub    = InterPolSub*v1;
%         v2Sub    = InterPolSub*v2;
%         vSub = [v1Sub;v2Sub];
% 
%         optsPlot.faceColour={'b'};
%         optsPlot.lineColour={'b'};
%         optsPlot.lineStyle=lineStyleDDFT{iDDFT};
% 
%         %plotRhoVdistDDFT2D(rhotSub,vSub,subArea,optsPlot,handlesRP);

        subMasses(iPlot,iDDFT) = subArea.Int*rhotSub;

    end
    
    
    lineColour = lineColourDDFT{iDDFT};
    lineColour = lineColour{1};
    
    plot(hMa,plotTimes,100*subMasses(:,iDDFT)/subMasses(1,iDDFT),'Color',lineColour)
    %plot(hMa,plotTimes,100*subMasses(:,iDDFT)/subMasses(end,iDDFT),'Color',lineColour)
    
    hold(hMa,'on');
    
end
    
yLabel = ['Percentage original mass within ' num2str(subArea.y2Max) ' $\sigma$ of wall'];

fixPlot(hMa,plotTimes(1),plotTimes(end),60,100,'Time',yLabel,[],'off',[]);
set(get(hMa,'YLabel'),'Rotation',90)

addOneLeg(optsPlot,hMf);

%write the figure files
save2pdf(outputFile,hMf);
close(hMf);

outputFiles = outputFile;


    %----------------------------------------------------------------------
    % Set axes, legend, add time
    %----------------------------------------------------------------------
    
    
    
    
%     optsPlot.time=[];
%     
%     optsPlot.xLab='Time';
%     optsPlot.yLab='y';
%     
%     optsPlot.xMin=optsPlot.rMin(1);
%     optsPlot.xMax=optsPlot.rMax(1);
%     optsPlot.yMin=optsPlot.rMin(2);
%     optsPlot.yMax=optsPlot.rMax(2);
% 
%     optsPlot.zMin=optsPlot.RMin;
%     optsPlot.zMax=optsPlot.RMax;
% 
%     optsPlot.zLab='Density';
% 
%     %optsPlot.legText=optsPlot.legTextR{iSpecies};
%     if(strcmp(optsPlot.plotType,'surf'))
%         fixPlot2Dsurf(hRa,optsPlot);
%         optsPlot.time=[];
%         fixPlot2Dcontour(hPa,optsPlot);
%     else
%         fixPlot2Dcontour(hRa,optsPlot);
%     end
%         
%     if(doP)  
% 
%         optsPlot.xMin=optsPlot.pMin(1);
%         optsPlot.xMax=optsPlot.pMax(1);
%         optsPlot.yMin=optsPlot.pMin(2);
%         optsPlot.yMax=optsPlot.pMax(2);            
% 
%         %optsPlot.legText=optsPlot.legTextP{iSpecies};
%         optsPlot.zLab='Momentum';
%         optsPlot.zMin=optsPlot.PMin(1);
%         optsPlot.zMax=optsPlot.PMax(1);
% 
%         fixPlot2D(hPa,optsPlot);
%     end
        
    
%     % write the figure files
%     save2pdf(outputFile,hRPf,100,true);
%     %close(hRPf);
%     
%     outputFiles = cat(2,outputFiles,outputFile);

fprintf(1,'Finished\n');

end % for iPlot


