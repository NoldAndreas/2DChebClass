function PlotThetaVsR()

    if(exist('D:\','dir'))
        dir = 'D://SyncGit/Projects/Data/ExperimentalContactAngle/';
    elseif(exist('/Users/NoldAndreas/','dir'))
        dir = '/Users/NoldAndreas/Documents/SyncGit/Projects/Data/ExperimentalContactAngle/';
    end

    nData = 1;
    data{1} = LoadDataFluidData('RameGaroff1996/Ca_0_005','s','r',struct('rad','deg','Ca','V')); 
	data{2} = LoadDataFluidData('RameGaroff1996/Ca_0_1','s','r',struct('rad','deg','Ca','V')); 

    PlotData(data{2});
    
    function PlotData(dat)
        figure('color','white','Position',[0 0 800 800]);
        rT = dat.r/dat.a;
        
        plot(rT,dat.theta,dat.symbol,'MarkerSize',10,'MarkerFaceColor',dat.color,'MarkerEdgeColor',dat.color); hold on;
        
        rTAna = min(rT) + (max(rT)-min(rT))*(0:0.01:1)';
        thetaAna = 180/pi*GHR_Inv(GHR_lambdaEta0(dat.thetaApp*pi/180)+ dat.Ca*log(rTAna),0);
        plot(rTAna,thetaAna,['-k'],'linewidth',1.5);
       % xlim([1e-5 100]);        ylim([0 180]);
        xlabel('$r/a$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta [^\circ]$','Interpreter','Latex','fontsize',20);

        set(gca,'fontsize',20);
    end
    
    function dataL = LoadDataFluidData(name,symbol,color,opts)
        fid = fopen([dir,name,'.txt']);
        y = textscan(fid,'%[^\n]',1,'headerlines',7); %[T, rhoG, rhoL]        
        x = textscan(fid,'%f %f'); %[T, rhoG, rhoL]
        fclose(fid); 
        
        
        dataL.Ca       = x{1}(1);
        dataL.a        = x{1}(2); % capillary length
        dataL.thetaApp = x{1}(3);
        dataL.thetaEq  = x{1}(4);
        
        dataL.r     = x{1}(4:end);
        
        if(strcmp(opts.rad,'rad'))
            dataL.theta = x{2}(4:end)*180/pi;            
        elseif(strcmp(opts.rad,'deg'))
            dataL.theta = x{2}(4:end);
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
end