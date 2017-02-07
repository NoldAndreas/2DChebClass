function str = getPDFMovieFile(filename,i,path)

    global dirData

    if(ispc())
        addPath = '\MoviesAuxPdfs\';
    else
        addPath = '/MoviesAuxPdfs/';
    end
    
    if(nargin == 3 && path == true)
        str = ([dirData,addPath]);
    else    
        %time = clock; % Gets the current time as a 6 element vector
        dirname = [dirData,addPath];        
        str = ([dirname,filename,'_',num2str(i),'.pdf']);
        if(~exist(dirname,'dir'))
            mkdir(dirname);
        end
    end
        
end