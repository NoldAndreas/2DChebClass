fileList = dir('*.fig');

nFiles = length(fileList);

for iFile = 1:nFiles
    
    fileName = fileList(iFile).name;
    open(fileName);
    pdfName = [fileName(1:end-3) 'pdf'];
    save2pdf(pdfName,gcf);
    
end