function combinePdfs(outputFile,pdfFileNames,quiet)
%combinePdfs(outputFile,pdfFileNames,quiet)
%  combines pdf files listed in pdfFileNames into outputFile
%
% INPUTS:
%  outputFile       --  a string for the output pdf file
%  pdfFileNames     --  a string list of the pdf file names
%  quiet            --  true/false whether to suppress output of gs

% if(quiet)
%     gsCmd= ['gs -dNOPAUSE -sDEVICE=pdfwrite ' ...
%         '-sOUTPUTFILE=' outputFile ' -dBATCH -dQUIET ' pdfFileNames];
% else
%     gsCmd= ['gs -dNOPAUSE -sDEVICE=pdfwrite ' ...
%         '-sOUTPUTFILE=' outputFile ' -dBATCH ' pdfFileNames];
% end

% this seems to work much better than gs
gsCmd=['pdftk ' pdfFileNames ' cat output ' outputFile];

system(gsCmd);
