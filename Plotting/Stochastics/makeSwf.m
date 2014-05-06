function makeSwf(swfFile,pdfFileNames,optsPlot)
%makeSwf(swfFile,pdfFileNames,optsPlot)
% combine the pdf files listed in pdfFileNames into swfFile
%
% INPUTS:
%  swfFile       --  a string for the swf file
%  pdfFileNames  --  a string list of the pdf file names
%  optsPlot      --  a size [1 1] structure containing
%             pdfDir          (string for directory in which to save pdfs)
%             movieFile       (string for output file name (NO EXTENSION) )
%             fps             (integer frames per second)
%             dpi             (integer resolution)
%             bitmap          (true/false whether to flatten pdfs -- needed
%                              for large animations)
%             quiet           (true/false whether to suppress output of gs
%                              and pdf2swf)


% get options
pdfDir=optsPlot.pdfDir;
fps=optsPlot.fps;
dpi=optsPlot.dpi;
bitmap=optsPlot.bitmap;
quiet=optsPlot.quiet;

fprintf(1,'Making pdf file ... ');

% pdf file in which to collate all time steps
pdfFile=[pdfDir 'All.pdf'];

% combine the pdf files into one
combinePdfs(pdfFile,pdfFileNames,quiet)

fprintf(1,'Finished\n');

% make swf
fprintf(1,'Making swf ... ');

% if we want to bitmap the pdf files (for flattening)
if(bitmap)
    swfCmd=['pdf2swf ' pdfFile ' -o ' swfFile  ' -s poly2bitmap -s multiply=1'...
            ' -s framerate=' num2str(fps) ' -s zoom=' num2str(dpi)];
else
    swfCmd=['pdf2swf ' pdfFile ' -o ' swfFile  ...
            ' -s framerate=' num2str(fps) ' -s zoom=' num2str(dpi)];
end

% if we want to suppress output of pdf2swf
if(quiet)
    swfCmd=[swfCmd ' -q'];
end
system(swfCmd);


fprintf(1,'Finished\n');