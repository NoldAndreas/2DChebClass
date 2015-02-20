function PlotThetaVsR()

    

    if(exist('D:\','dir'))
        dir = 'D://SyncGit/Projects/Data/ExperimentalContactAngle/';
    elseif(exist('/Users/NoldAndreas/','dir'))
        dir = '/Users/NoldAndreas/Documents/SyncGit/Projects/Data/ExperimentalContactAngle/';
    end
    ChangeDirData(dir,'Org');
    
 %   LoadDataFluidData('\ChenRameGaroff/Fig4a','s','r',struct('rad','deg','Ca','V')); 
 %   LoadDataFluidData('\ChenRameGaroff/Fig6a','s','r',struct('rad','deg','Ca','V')); 
    data = LoadDataFluidData('RameGaroff1996/Ca_0_005','s','b',struct('rad','deg','Ca','V')); 
    PlotData(data,'RameGaroff1996/RameGaroff1996_Ca_0_005');
	data = LoadDataFluidData('RameGaroff1996/Ca_0_1','s','b',struct('rad','deg','Ca','V')); 
    PlotData(data,'RameGaroff1996/RameGaroff1996_Ca_0_1');
    

    %PlotData(data{1});
    
    function PlotData(dat,filename)
        figure('color','white','Position',[0 0 800 800]);
        rT = dat.r/dat.a;
        
        plot(rT,dat.theta,dat.symbol,'MarkerSize',6,'MarkerFaceColor',dat.color,'MarkerEdgeColor',dat.color); hold on;
        
        rTAna   = min(rT) + (max(rT)-min(rT))*(0:0.01:1)';        
        f_0     = f_0_HuhScrivenCapillary(rTAna,dat.thetaAppRad,pi/2,dat.R_T/dat.a);%pi/2
        f_0DR   = f_0_DussanRameGaroff(rTAna,dat.thetaAppRad,pi/2,dat.R_T/dat.a);%pi/2
        thetaCox = GHR_Inv(GHR_lambdaEta0(dat.thetaAppRad)+ dat.Ca*log(rTAna),0);
        
        thetaAna = thetaCox + f_0 - dat.thetaAppRad;
        
        plot(rTAna,180/pi*f_0,['--k'],'linewidth',1.5);
        %plot(rTAna*dat.a,180/pi*f_0DR,['--b'],'linewidth',1.5);        
        plot(rTAna,180/pi*thetaAna,['-k'],'linewidth',1.5);        
       % xlim([1e-5 100]);        ylim([0 180]);
       
       
        xlim([min(rT),max(rT)]);
        xlabel('$r/a$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta [^\circ]$','Interpreter','Latex','fontsize',20);
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        
        SaveFigure([filename,'_fig'],dat);
    end
    
    function dataL = LoadDataFluidData(filename,symbol,color,opts)
        fid = fopen([dir,filename,'.txt']);
        y = textscan(fid,'%[^\n]',1,'headerlines',8); %[T, rhoG, rhoL]        
        x = textscan(fid,'%f %f'); %[T, rhoG, rhoL]
        fclose(fid); 
        
        
        dataL.Ca       = x{1}(1);
        dataL.a        = x{1}(2); % capillary length
        dataL.thetaAppDeg = x{1}(3);
        dataL.thetaAppRad = dataL.thetaAppDeg*pi/180;
        dataL.thetaEq  = x{1}(4);
        dataL.R_T      = x{1}(5);
        
        dataL.r     = x{1}(6:end);
        
        if(strcmp(opts.rad,'rad'))
            dataL.theta = x{2}(6:end)*180/pi;            
        elseif(strcmp(opts.rad,'deg'))
            dataL.theta = x{2}(6:end);
        end
        
        dataL.legend  = char(y{1});        
        dataL.color   = color;
        dataL.symbol  = symbol;
        
        if(isempty(color))
            dataL.lw = 2;            
        else
            dataL.MarkerSize = 10;            
        end   
                 
    end
    function th = f_0_HuhScrivenCapillary(r,thOut,alpha,R_T)
        phi = 1/4*(alpha - thOut);
        phi_r = atan(besselk(1,r+R_T)/besselk(1,R_T)*tan(phi));  
        th  = alpha - 4*phi_r;   
        th  = real(th);
    end
end