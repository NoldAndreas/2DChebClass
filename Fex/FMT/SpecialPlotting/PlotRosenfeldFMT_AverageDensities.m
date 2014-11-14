function PlotRosenfeldFMT_AverageDensities(FMTShape,FMTMatrices,rho)

    figure('name','Average densities'); 
    set(gcf, 'Position', [0 0 1500 1000]);    set(gcf,'Color','white');
    nStruct = FMTMatrices.AD;    
    fields  = fieldnames(nStruct);
    
    noRows  = ceil(sqrt(numel(fields)));
    noCols  = ceil(numel(fields)/noRows);

    for i=1:numel(fields)
      %fields(i)
      %nStruct.(fields(i))      
        subplot(noRows,noCols,i);
        FMTShape.AD.plot(nStruct.(fields{i})*rho);    
        title(fields(i));    
    end

%     subplot(2,2,1);
%     FMTShape.PlotFull(FMTMatrices.AD.n2*rho,1);    
%     title('n_2');
%     
% 	subplot(2,2,2);
%     FMTShape.PlotFull(FMTMatrices.AD.n1*rho,1);
%     title('n_1');
%     
%     subplot(2,2,3);
%     FMTShape.PlotFull(FMTMatrices.AD.n1_v_1*rho,1);    
%     title('n1_v_1');    
%     
%     subplot(2,2,4);
%     FMTShape.PlotFull(FMTMatrices.AD.n1_v_2*rho,1);       
%     title('n1_v_2');                

end