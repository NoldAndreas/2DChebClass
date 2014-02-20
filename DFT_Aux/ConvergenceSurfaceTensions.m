function ConvergenceSurfaceTensions(config)

    config.optsNum.maxComp_y2         = -1;
  %  config.optsNum.PhysArea.alpha_deg = 90;
    
    ST50 = ComputeNormal(50);
    ST60 = ComputeNormal(60);
    ST70 = ComputeNormal(70);
    ST80 = ComputeNormal(80);
    ST90 = ComputeNormal(90);
    ST100 = ComputeNormal(100);  
    
    
    %********************************************
    config.optsNum.PhysArea.N2bound = 3;
     
    STP30 = ComputeParallel(30);
    STP40 = ComputeParallel(40);    
	STP50 = ComputeParallel(50);
    STP60 = ComputeParallel(60);
    STP70 = ComputeParallel(70);
           
    format long
    wl_ref    = ST80.om_wallLiq;
    wg_ref    = ST80.om_wallGas;
    lg_ref    = STP40.om_LiqGas;
    theta_ref = ComputeContactAngle(wg_ref,wl_ref,lg_ref);
    
    wg_dat = [ST50.om_wallGas,ST60.om_wallGas,ST70.om_wallGas,ST80.om_wallGas,ST90.om_wallGas];%,ST100.om_wallGas,ST110.om_wallGas];
    wl_dat = [ST50.om_wallLiq,ST60.om_wallLiq,ST70.om_wallLiq,ST80.om_wallLiq,ST90.om_wallLiq];%,ST100.om_wallLiq,ST110.om_wallLiq];
    lg_dat = [STP30.om_LiqGas,STP40.om_LiqGas,STP50.om_LiqGas];
    
    %Max Deviation:
    %cos(theta)-cos(theta_ref) = -sin(theta_ref)*(theta-theta_ref) + HOT
    
    %Max Deviation due to Wall-Gas Inacurracy
    wg_theta_err = 1/sin(theta_ref)*max(abs(wg_dat - wg_ref))/lg_ref;
    wl_theta_err = 1/sin(theta_ref)*max(abs(wl_dat - wl_ref))/lg_ref;
    lg_theta_err = 1/sin(theta_ref)*max(abs(1./lg_dat - 1./lg_ref))*abs(wg_ref-wl_ref);
    
    fprintf('************************************************\n');
    fprintf(['Max Error due to wg-inacurracies: ',num2str(wg_theta_err*180/pi),'[deg]\n']);
    fprintf(['Max Error due to wl-inacurracies: ',num2str(wl_theta_err*180/pi),'[deg]\n']);
    fprintf(['Max Error due to lg-inacurracies: ',num2str(lg_theta_err*180/pi),'[deg]\n']);
    fprintf('************************************************\n');
    
    %CA comparison
     alpha1 = ComputeContactAngle(ST50.om_wallGas,...
                                   ST50.om_wallLiq,...
                                   STP30.om_LiqGas)
%                               
     alpha2 = ComputeContactAngle(ST70.om_wallGas,...
                                   ST70.om_wallLiq,...
                                   STP50.om_LiqGas)                              
%                               
%     alpha100 = ComputeContactAngle(ST100.om_wallGas,...
%                                    ST100.om_wallLiq,...
%                                    STP100.om_LiqGas);
                               
    function ST = ComputeNormal(N2)
        config.optsNum.PhysArea.N = [1;N2];
        CL = ContactLine(config); %CL with normal high resolution
        CL.Preprocess();     
        %CL.ComputeST(true);
        ST = CL.ST_1D;
    end

    function ST = ComputeParallel(N1)
        config.optsNum.PhysArea.N = [N1;3];
        CL = ContactLine(config); 
        CL.Preprocess();        
       % CL.ComputeST(true);
        ST = CL.ST_1D;                
    end
       

end