function optsNumDDFT=getOptsDDFT(optsStruct)
% optsStoc=getOptsDDFT(optsStruct,fileStruct)
%  Collates DDFT options into the output structure
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   doDDFT:        t/f do DDFT calculation {nDDFT,1}
%                   loadDDFT:      t/f load DDFT calculation {nDDFT,1}
%                   saveDDFT:      t/f save DDFT calculation {nDDFT,1}
%                   DDFTtype:      'r'/'rv' without/with inertia {'x1',...,'xnDDFT')
%                   DDFTName:      name of DDFT calculation {'name1',...,'namenDDFT') 
%                   DDFTDir:       directory to add to path (string or 
%                                    {string 1, ..., string nDDFT}
%                   DDFTCode:      names of DDFT code (string or 
%                                    {string 1, ..., string nDDFT}
%                   DDFTParams     struct of size [1 nDDFT] containing DDFT
%                                      parameters
%
%   fileStruct  -- struct of size [1 1] containing
%                   saveFileDDFT:  string save file for DDFT data
%
% OUTPUTS:
%   optsStoc  -- struct of size [1 1] containing all the above data
if( optsStruct.anyDDFT )   
    optsNumDDFT=struct('plotTimes',optsStruct.plotTimes, ...
                     'saveDDFT',optsStruct.saveDDFT, ...
                     'loadDDFT',optsStruct.loadDDFT, ...
                     'doDDFT',optsStruct.doDDFT, ...
                     'type',optsStruct.DDFTType,...
                     'DDFTName',optsStruct.DDFTName, ...
                     'DDFTCode',optsStruct.DDFTCode);
                     %'DDFTDir',optsStruct.DDFTDir, ...
                     %'saveFileDDFT',fileStruct.saveFileDDFT); 
                    
    % merge in the user-specified DDFT parameters
    DDFTParams=optsStruct.DDFTParams;
    
    optsNumDDFT=mergeStruct(optsNumDDFT,DDFTParams); 

    % cut down to the calculations we want to do
    optsNumDDFT=optsNumDDFT(cell2mat(optsStruct.doDDFT));
    
else
    optsNumDDFT=[];
end

end