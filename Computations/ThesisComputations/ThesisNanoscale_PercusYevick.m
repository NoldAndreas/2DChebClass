function ThesisNanoscale_PercusYevick

    AddPaths('ThesisNanoscale');
    filename = 'D:\Data\PercusYevick\PY05.txt';
    
	fid = fopen(filename);
	x   = textscan(fid,'%f %f'); %[T, rhoG, rhoL],'headerlines',3
	fclose(fid); 
	
    figure('Position',[0 0 300 250],'Color','white');
    
    plot(x{1},x{2},'k','MarkerFaceColor','k'); hold on;	
    xlim([0 4]);
    xlabel('$r$','Interpreter','Latex');
    ylabel('$\radialDistributionFunction(r)$','Interpreter','Latex');
    SaveFigure('PercusYevick_RadialDistributionFunciton');
    
	%plot(x{3},x{1},['k',sym],'MarkerFaceColor','k'); hold on;	       

end