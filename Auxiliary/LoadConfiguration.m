function [configuration,configName] = LoadConfiguration()
    global dirData
    
    [configIn,DataFolder] = uigetfile([dirData filesep 'Configurations' filesep '*.mat'],['Select Config File']);
    load([DataFolder,configIn]);
    disp(['Loading configuration ',[DataFolder,configIn],' ...']);
    configName = configIn(1:end-4);        
end