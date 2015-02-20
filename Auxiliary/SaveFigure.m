function fullName = SaveFigure(filename,opts)

    global dirData
    
    k = strfind(filename,'.');
    while(~isempty(k))
        filename(k) = '_';
        k = strfind(filename,'.');
    end
    
    [s,branch]= system('C:\git rev-parse --abbrev-ref HEAD');
	
    fullName = [dirData filesep filename];
    
    DataFolder = fileparts(fullName);
	if(~exist(DataFolder,'dir'))            
        disp('Folder not found. Creating new path..');            
        mkdir(DataFolder);
	end
    
    print2eps(fullName,gcf);
	saveas(gcf,[fullName '.fig']);        
    
    disp(['Figures saved in ',fullName '.fig/eps']);
        
    if((nargin >= 2) && ~isempty(opts)) 
        Struct2File([fullName, '.txt'],opts,...
                    ['Figure ',filename,'.fig/eps saved at: ',datestr(now),' with git branch ',branch(1:end-1)]);
    end
            
	disp(['Options saved in ',fullName '.tex']);

end