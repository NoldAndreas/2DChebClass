function epw = FindEpwFromContactAngle(config,thY)
    
    
    %Get liquid-vapour surface tension
    config.optsPhys.V1.epsilon_w = 0;
    config.optsNum.maxComp_y2    = -1;        
    
    confP = config;
    confP.optsNum.PhysArea.N2bound = 3;
    confP.optsNum.PhysArea.N = [60;3];    
    
    CLP = ContactLineHS(confP);
    CLP.Preprocess();     close all;
    omLG = CLP.ST_1D.om_LiqGas;            
    
    %**********************
    config.optsNum.PhysArea.N = [1;100];
	CLN = ContactLineHS(config);
    CLN.Preprocess();     
    close all;        
    
    epw     = zeros(size(thY));
    options = optimoptions('fsolve','TolX',0.0001);%'Display','off');
    
    for i = 1:length(thY)
        thY_rad_i = thY(i)*pi/180;        
        epw(i)  = fsolve(@FindEpw,0.9,options);            
    end
        
    
    function z = FindEpw(epw)                
        
        CLN.optsPhys.V1.epsilon_w = epw;      
        CLN.Compute1D('WL'); close all;
        CLN.Compute1D('WG'); close all;
            
        omWG = CLN.ST_1D.om_wallGas;
        omWL = CLN.ST_1D.om_wallLiq;
                        
        z = (omWG - omWL)/omLG - cos(thY_rad_i);

    end
end