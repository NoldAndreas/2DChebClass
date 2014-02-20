function str = getPDFMovieFile(filename,i,path)

    if(ispc())
        addPath = 'Movies\PDFs\';
    else
        addPath = 'Movies/PDFs/';
    end
    
    if(nargin == 3 && path == true)
        str = ([getResultsPath(),addPath]);
    else    
        %time = clock; % Gets the current time as a 6 element vector
        str = ([getResultsPath(),addPath,filename,'_',num2str(i),'.pdf']);
    end
        
end