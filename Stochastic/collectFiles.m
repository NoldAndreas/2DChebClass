function collectFiles(inputFile,outputDir)
    if(~exist(inputFile,'file'))
        disp('File does not exist');
        return;
    end
    
    if(~exist(outputDir,'dir'))
        mkdir(outputDir);
    end
    
    [~,inputFileName,inputFileExt] = fileparts(inputFile);
    inputFileName = [inputFileName inputFileExt];
    copyfile(inputFile,[outputDir filesep inputFileName]);
    
    fid = fopen(inputFile,'r');
    
    tline = fgets(fid);
    while ischar(tline)
        mFile = [tline(2:end-1) '.mat'];
        if(exist(mFile,'file'))
            disp(mFile)
            [~,mFileName,~] = fileparts(mFile);
            copyfile(mFile,[outputDir filesep mFileName '.mat'])
            tFile = [tline(2:end-1) '.txt'];
            copyfile(tFile,[outputDir filesep mFileName '.txt'])
        end
        tline = fgets(fid);
    end

    fclose(fid);
end