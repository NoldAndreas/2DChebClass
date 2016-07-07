fileList = dir('*.fig');

nFiles = length(fileList);
pdfFileNames = [];

fprintf(1,'Making pdf ... ');

for iFile = 1:nFiles
    
    fileName = fileList(iFile).name;
    pdfName = [fileName(1:end-3) 'pdf'];
    pdfFileNames = cat(2, pdfFileNames, [' ' pdfName]);
    if(~exist(pdfName,'file'))
        h = open(fileName);
        save2pdf(pdfName,gcf);
        close(h);
    end
    
end

fprintf(1,'Combining pdf ... ');
fullPdfFile=['movie.pdf'];

gsCmd= ['gs -dNOPAUSE -sDEVICE=pdfwrite ' ...
          '-sOUTPUTFILE=' fullPdfFile ' -dBATCH -dQUIET ' pdfFileNames];

system(gsCmd);
fprintf(1,'Finished\n');