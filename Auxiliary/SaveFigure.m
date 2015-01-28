function SaveFigure(filename,opts)

    global dirData
    
    k = strfind(filename,'.');
    while(~isempty(k))
        filename(k) = '_';
        k = strfind(filename,'.');
    end
    
    [s,branch]= system('C:\git rev-parse --abbrev-ref HEAD');
	
    print2eps([dirData filesep filename],gcf);
	saveas(gcf,[dirData filesep filename '.fig']);        
    
    disp(['Figures saved in ',dirData filesep filename '.fig/eps']);
        
    if((nargin >= 2) && ~isempty(opts)) 
        Struct2File([dirData filesep filename, '.txt'],opts,...
                    ['Figure ',filename,'.fig/eps saved at: ',datestr(now),' with git branch ',branch(1:end-1)]);
    end
            
	disp(['Options saved in ',dirData filesep filename '.tex']);

end