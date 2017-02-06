function plotInitialFinalParticles(meanStruct,optsPlot,gifFile)
%plotInitialFinalParticles(meanStruct,optsPlot,gifFile)
% plots initial and final mean positions of particles given a mean
% structure made by getMeansParticles
%
% INPUTS:
%  meanStruct   --  structure of size [1 nStoc] containing
%                    xMean         (matrix of size 3*nParticles x nTimes)
%                    pMean         (matrix of size 3*nParticles x nTimes)
%
%  optsPlot     --  structure of size [1 1] containing
%                    sigma         (particle diameter)
%                    colours       (colours for each plot, of the form 
%                                    {{'r',...,'b'}})
%                    shift         (3xnPlots matrix to shift the plots by)
%                    renormalize   (3x1 vector of true/false whether to 
%                                    renormalize coordinate)
%                    relabel       (3x1 vector of true/false whether to 
%                                    relabel axes)
%                    lims          (cell of three 1x2 vectors of limits,
%                                 {{[xmin,xmax],[ymin,ymax],[zmin,zmax]}} )
%                    ticks         (cell of tick mark locations given 
%                                     by row vectors, e.g.{{[0 2],[],[]}} )
%                    labs          (cell of list of tick mark labels given 
%                                  as column vectors,
%                                  e.g. {{{'HI off'; 'HI on'},{''},{''}}} )
%                    followType    ('max'/'min', whether to follow the max or
%                                 min of a coordinate when using renormalize)
%
%  gifFile    --     1x3 cell of strings for output files (only last two
%                     are used)

fprintf(1,'Making particle initial and final plots ... ');

% number of plots to do at each time
nStoc=size(meanStruct,2);

% number of coordinates = 3*nParticles
meanTemp=meanStruct(1).xMean;
nCoords=size(meanTemp,1);

% initialize matrices
xInitial=zeros(nCoords,nStoc);
xFinal=xInitial;
pInitial=xInitial;
pFinal=xInitial;

% get initial and final data for each plot
for iStoc=1:nStoc
    xInitial(:,iStoc)=meanStruct(iStoc).xMean(:,1);
    xFinal(:,iStoc)=meanStruct(iStoc).xMean(:,end);
    
    pInitial(:,iStoc)=meanStruct(iStoc).pMean(:,1);
    pFinal(:,iStoc)=meanStruct(iStoc).pMean(:,end);
end

pScale=optsPlot.pScale;
pInitial=pInitial/pScale;
pFinal=pFinal/pScale;

% set up figure and axes
fullscreen = get(0,'ScreenSize');
hIFf=figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
%hIFf=figure;
hIa=subplot(1,2,1);
hFa=subplot(1,2,2);

% collect handles
handlesI=struct('hf',hIFf,'ha',hIa);
handlesF=struct('hf',hIFf,'ha',hFa);

% set background to be white
set(hIFf, 'Color', 'w');

% do initial plot
optsPlot.title='Initial';
plotParticles3D(xInitial,pInitial,optsPlot,handlesI);
% do final plot
optsPlot.title='Final';
plotParticles3D(xFinal,pFinal,optsPlot,handlesF);

fprintf(1,'Saving ...\n');
%save2pdf(pdfFile,hIFf,optsPlot.dpi);
%save2tiff('temp.tiff',hIFf,optsPlot.dpi);
f = getframe(hIFf);
[im,map] = rgb2ind(f.cdata,256,'nodither');
imwrite(im,map,gifFile);

fprintf(1,'Saved\n');

%close(hIFf);

fprintf(1,'Finished\n');