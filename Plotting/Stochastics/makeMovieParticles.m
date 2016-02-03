function makeMovieParticles(meanStruct,optsPlot)
%makeMovieParticles(meanStruct,optsPlot)
% makes movie of mean positions of particles given a mean structure made by
% getMeansParticles
%
% INPUTS:
%  meanStruct   --  structure of size [1 nStoc] containing
%                    xMean         (matrix of size 3*nParticles x nTimes)
%                    pMean         (matrix of size 3*nParticles x nTimes)
%
%  optsPlot     --  structure of size [1 1] containing
%             sigma         (particle diameter)
%             colours       (colours for each plot, of the form 
%                            {{'r',...,'b'}})
%             shift         (3xnPlots matrix to shift the plots by)
%             renormalize   (3x1 vector of true/false whether to 
%                            renormalize coordinate)
%             relabel       (3x1 vector of true/false whether to 
%                            relabel axes)
%             lims          (cell of three 1x2 vectors of limits,
%                         {{[xmin,xmax],[ymin,ymax],[zmin,zmax]}} )
%             ticks         (cell of tick mark locations given 
%                             by row vectors, e.g.{{[0 2],[],[]}} )
%             labs          (cell of list of tick mark labels given 
%                          as column vectors,
%                          e.g. {{{'HI off'; 'HI on'},{''},{''}}} )
%             followType    ('max'/'min', whether to follow the max or
%                         min of a coordinate when using renormalize)
%             plotTimes     (vector  of plotting times)
%             doMovieGif      (true/false whether to make gif)
%             doPdfs          (true/false whether to save pdfs for swf)
%             doMovieSwf      (true/false whether to make swf)
%             pdfDir          (string for directory in which to save pdfs)
%             movieFile       (string for output file name (NO EXTENSION) )
%             fps             (integer frames per second)
%             dpi             (integer resolution)
%             bitmap          (true/false whether to flatten pdfs -- needed
%                              for large animations)
%             quiet           (true/false whether to suppress output of gs
%                              and pdf2swf)

fprintf(1,'Making particle movie plots ... ');

optsPlot=optsPlot(1);

% set up figure and axes
fullscreen = get(0,'ScreenSize');
%hRf=figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
hRf=figure('Position',[0 0 600 600]);
%hRf=figure;
hRa=axes;
% collect handles
handles=struct('hf',hRf,'ha',hRa);

% stops flashing
set(hRa,'nextplot','replacechildren');

% set background to be white
set(hRf, 'Color', 'w');

% get number of data sets at each time
nStoc=size(meanStruct,2);
% and number of coordinates = 3*nParticles
meanTemp=meanStruct(1).xMean;
nCoords=size(meanTemp,1);

% get plot times
plotTimes=optsPlot.plotTimes;
nPlots=length(plotTimes);

% initialize matrices
x=zeros(nCoords,nStoc);
p=x;

doGif=optsPlot.doMovieGifP;
doSwf=optsPlot.doMovieSwfP;
doPdfs=optsPlot.doPdfsP;

% get movie files and options
pdfDir=optsPlot.pdfDir;
movieFile=optsPlot.movieFile;
fps=optsPlot.fps;
dpi=optsPlot.dpi;

% set up pdf file names
if(doPdfs || doSwf)
    if(~exist(pdfDir,'dir'))
        mkdir(pdfDir);
    end       
    pdfFileNames=[];
    nDigits=ceil(log10(nPlots));
    nd=['%0' num2str(nDigits) 'd'];
end

for iPlot=1:nPlots
    % clear axis
    cla(hRa);
    
    % get positions for current time
    for iStoc=1:nStoc
        x(:,iStoc)=meanStruct(iStoc).xMean(:,iPlot);
        p(:,iStoc)=meanStruct(iStoc).pMean(:,iPlot);
    end

    % scale the momentum for arrows
    p=p/optsPlot.pScale;
    
    % set title
    %optsPlot.title=['t=' num2str(plotTimes(iPlot))];
    optsPlot.title=[];
    
    % plot particles
    plotParticles3D(x,p,optsPlot,handles)
     
    if(doGif)
        % plot the current frame
        f = getframe(hRf);

        % get image data
        if(iPlot==1)
            % think this makes sure that the color map is correct,
            % otherwise it is black and white
            [im,map] = rgb2ind(f.cdata,256,'nodither');

            % preallocate something in the final frame which preallocates the
            % entire matrix, giving a significant increase in performance
            im(1,1,1,nPlots) = 0;
        else
            im(:,:,1,iPlot) = rgb2ind(f.cdata,map,'nodither');
        end
    end
    
    if(doPdfs || doSwf)
        outputFile=[pdfDir num2str(iPlot,nd) '.pdf'];
    
        pdfFileNames = cat(2, pdfFileNames, [' ' outputFile]);
    
        if(doPdfs)
            save2pdf(outputFile,hRf,dpi);
        end
              
    end
    
end

fprintf(1,'Finished\n');

% save the file

if(doGif)
    fprintf(1,'Making gif movie ... ');
    % write the movie file
    gifFile=[movieFile '.gif'];
    delayTime=1/fps;
    %imwrite(im,map,gifFile,'DelayTime',delayTime,'LoopCount',inf)
    imwrite(im,map,gifFile,'DelayTime',delayTime,'LoopCount',0)
    fprintf(1,'Finished\n');
end

if(doSwf)
    swfFile=[movieFile '.swf'];
    fprintf(1,'Making swf movie ... ');
    swfCmd=['gif2swf ' gifFile ' -o ' swfFile];
    system(swfCmd);
    fprintf(1,'Finished\n');
end

if(doPdfs)
    fprintf(1,'Combining pdf ... ');
    fullPdfFile=[movieFile '.pdf'];

    gsCmd= ['gs -dNOPAUSE -sDEVICE=pdfwrite ' ...
              '-sOUTPUTFILE=' fullPdfFile ' -dBATCH -dQUIET ' pdfFileNames];

    system(gsCmd);
    fprintf(1,'Finished\n');
end

% if(doSwf)
%     % name for swf file
%     swfFile=[movieFile '.swf'];
%     
%     % make swf
%     makeSwf(swfFile,pdfFileNames,optsPlot)
% end

close(hRf);