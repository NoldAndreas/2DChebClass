function addPaths()
    % addPaths()
    %   Adds the necessary paths in the folder containing dynamicsComparison
    %
    % INPUTS:
    %   <none>
    %
    % OUTPUTS:
    %   <none>

    path('Inputs',path)           % input files
    path('PreProcessor',path)       % preprocessor files
    path('Potentials',path)       % potential files
    path('Stochastic',path)       % stochastic calculations
    path(['Stochastic' filesep 'InitialGuess'],path)     % initial guess files
    path(['Stochastic' filesep 'HI'],path)               % hydrodynamic interaction files
    path('Plotting',path)         % plotting

    global dirData

    if(exist('/home/bgoddard/','dir'))
        dirData = '/home/bgoddard/work/MATLAB/Fluids/NComponent2D/Data/2DChebData';
    elseif(exist('/Users/Ben/','dir'))
        dirData = '/Users/Ben/work/MATLAB/Fluids/NComponent2D/Data/2DChebData';
    else
        disp('Unknown computer; using current directory for save data');
    end

end