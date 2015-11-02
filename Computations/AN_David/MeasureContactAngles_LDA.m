function [optsNum,optsPhys] = MeasureContactAngles_LDA()
    
    close all;
    clear all;
    AddPaths(); 

    %Numerical Parameters    
    Phys_Area = struct('shape','HalfSpace_FMT','N',[1,20],...
                       'L1',4,'L2',5,'L2_AD',2.,...
                       'y2wall',0.,...
                       'N2bound',24,'h',1);
    
    Plot_Area = struct('y1Min',-10,'y1Max',10,'N1',100,'N2',100,...
                       'y2Min',0.5,'y2Max',10);
        
    Sub_Area = struct('shape','Box','y1Min',-1,'y1Max',1,'N',[20,20],...
                      'y2Min',0,'y2Max',2);
                                                        
	V2Num   = struct('Fex','SplitDisk','L',2,'L2',1.,'N',[20,20]);
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'plotTimes',0:0.1:6,...
                     'V2Num',V2Num);                                           
    
    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 
    
                 
    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','CarnahanStarling',...
                      'kBT',0.7,...                      
                      'Dmu',0.0,...
                      'sigmaS',1);
                  
    T   = 0.7;%0.6:0.05:0.85;
    epw = 1.5:0.2:3.2;
    c = {'k','r','b','m','g','y','k'};
    
    res = DataStorage([],@ComputeContactAngleVsEpsilon,v2struct(T,epw),[]);          
    
    f1 = figure('color','white','Position',[0 0 800 800]);
    for k = 1:length(T)
        plot(epw,res.theta(k,:)'*180/pi,c{k},'linewidth',1.5); hold on;
        plot(epw,res.theta_DWSI(k,:)'*180/pi,[c{k},'--'],'linewidth',1.5); hold on;
    end
    xlabel('$\varepsilon_w/\varepsilon$','Interpreter','Latex','fontsize',20);
    ylabel('$\theta [^\circ]$','Interpreter','Latex','fontsize',20);
    set(gca,'fontsize',20);
    set(gca,'linewidth',1.5);
    
    print2eps('ContactAngleVsEpsilon',f1);
    saveas(f1,'ContactAngleVsEpsilon.fig');
    
    f1 = figure('color','white','Position',[0 0 800 800]);
    for k = 1:length(T)
        plot(epw,res.om_wl(k,:)',[c{k}],'linewidth',1.5); hold on;
        plot(epw,res.om_wg(k,:)',[c{k},'--'],'linewidth',1.5); hold on;
    end
    xlabel('$\varepsilon_w/\varepsilon$','Interpreter','Latex','fontsize',20);
    ylabel('$\gamma$','Interpreter','Latex','fontsize',20);
    set(gca,'fontsize',20);
    set(gca,'linewidth',1.5);    
    
    print2eps('GammaWL',f1);
    saveas(f1,'GammaWL.fig');
    
    f1 = figure('color','white','Position',[0 0 800 800]);    
    plot(T,res.om_lg,'linewidth',1.5); hold on;            
    xlabel('$\varepsilon_w/\varepsilon$','Interpreter','Latex','fontsize',20);
    ylabel('$\gamma$','Interpreter','Latex','fontsize',20);
    set(gca,'fontsize',20);
    set(gca,'linewidth',1.5);
    
    print2eps('GammaLG',f1);
    saveas(f1,'GammaLG.fig');
    
    function res = ComputeContactAngleVsEpsilon(params,other)
        T   = params.T;
        epw = params.epw;
        
        theta     = zeros(length(T),length(epw));
        thetaDWSI = zeros(length(T),length(epw));
        
        om_wg = zeros(length(T),length(epw));
        om_wl = zeros(length(T),length(epw));        
        om_lg = zeros(length(T),1);    
                
        optsNum.PhysArea.N = [70,3];                
        EX                 = ContactLineHS(v2struct(optsPhys,optsNum));    
        EX.Preprocess();

        for j = 1:length(T)        
            EX.ResetTemperature(T(j));   
            om_lg(j)        = EX.Compute1D('LG');
        end

        optsNum.PhysArea.N = [1;70];
        EX     = ContactLineHS(v2struct(optsPhys,optsNum));
        EX.Preprocess();   

        for j = 1:length(T)

            EX.ResetTemperature(T(j));        

            for i = 1:length(epw)
                EX.optsPhys.V1.epsilon_w = epw(i);
                EX.Preprocess_ExternalPotential();
                om_wg(j,i) = EX.Compute1D('WG');
                om_wl(j,i) = EX.Compute1D('WL');                

                theta(j,i)     = ComputeContactAngle(om_wg(j,i),om_wl(j,i),om_lg(j));
                thetaDWSI(j,i) = DryWallSI(T(j),epw(i),...
                                           EX.optsPhys.rhoLiq_sat,...
                                           EX.optsPhys.rhoGas_sat);
                

                close all;                
                if(theta(j,i)==0)
                    break;
                end
            end
        end
        
        res.theta      = theta;
        res.theta_DWSI = thetaDWSI;
        res.om_lg      = om_lg;
        res.om_wg      = om_wg;
        res.om_wl      = om_wl;
        
    end
    
    
end        

function theta = DryWallSI(T,epw,nliq,ngas)
    theta = acos( 58/135*T*epw/(nliq-ngas) + (ngas^2-nliq^2)/(nliq-ngas)^2 );
end

