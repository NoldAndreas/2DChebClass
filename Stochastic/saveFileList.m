function saveFileList(inputFile,potNames,stocStruct,DDFTStruct)
    
    global dirData;
    inputFileCopy = [dirData filesep potNames filesep 'fileList_' getTimeStr() '.txt'];
    
    copyfile(which(inputFile),inputFileCopy);
    
    fid = fopen(inputFileCopy,'a');

    nStoc = length(stocStruct);
    nDDFT = length(DDFTStruct);
    
    divString = '%%--------------------------------------------------------------------------\n';
    
    if(nStoc > 0)
        fprintf(fid,'\n');
        fprintf(fid,divString);
        fprintf(fid,'%% Stochastic output files\n');
        fprintf(fid,divString);
        
        ICFilename = stocStruct(1).ICFilename;
        FCFilename = stocStruct(1).FCFilename;
        pFilename  = stocStruct(1).pFilename;
        
        fprintf(fid,'%% Sampling:\n');
        fprintf(fid,['%%' ICFilename '\n']);
        if(~isempty(FCFilename))
            fprintf(fid,['%%' FCFilename '\n']);
        else
            fprintf(fid,'%% no final sampling\n');
        end
        if(~isempty(pFilename))
            fprintf(fid,['%%' pFilename '\n']);
        else
            fprintf(fid,'%% no final sampling\n');
        end
        
        fprintf(fid,'%% Dynamics:\n');
    end

    for iStoc = 1:nStoc
        fprintf(fid,['%%' stocStruct(iStoc).Filename '\n']);
    end
    
    if(nDDFT > 0)
        fprintf(fid,'\n');
        fprintf(fid,divString);
        fprintf(fid,'%% DDFT output files\n');
        fprintf(fid,divString);
    end
    
    for iDDFT = 1:nDDFT
        fprintf(fid,['%%' DDFTStruct(iDDFT).Filename '\n']);
    end
    
    fclose(fid);
    
    
end
