function fullName = SaveFigure(filename,opts)

    global dirData
    
    k = strfind(filename,'.');
    while(~isempty(k))
        filename(k) = '_';
        k = strfind(filename,'.');
    end
    
    [s,branch]= system('C:\git rev-parse --abbrev-ref HEAD');
	
    if(isempty(fileparts(filename)));
        fullName = [dirData filesep filename];
    else
        fullName = filename;
    end
    
    DataFolder = fileparts(fullName);
	if(~exist(DataFolder,'dir'))            
        disp('Folder not found. Creating new path..');            
        mkdir(DataFolder);
    end
        
    print2eps(fullName,gcf);
	saveas(gcf,[fullName '.fig']);        
    matlab2tikz([fullName '.tex'],'showInfo',false);
    save2pdf([fullName '.pdf'],gcf);
    
    disp(['Figures saved in ',fullName '.fig/eps/tex/pdf']);
        
    if((nargin >= 2) && ~isempty(opts)) 
        opts_filename = [fullName, '_figData.txt'];
        Struct2File(opts_filename,opts,...
                    ['Figure ',filename,'.fig/eps saved at: ',datestr(now),' with git branch ',branch(1:end-1)]);
        disp(['Options saved in ',opts_filename]);
    end

end