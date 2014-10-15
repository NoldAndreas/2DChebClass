function [Pts,Diff,Int,Ind,Interp] = FDFD(Maps,N1,N2,interpPts1,interpPts2)
    % Variables in the computational domain [-1,1]: x_1, x_2
    % Variables in the physical domain: y_1, y_2    

    [x1,wInt1]  = FDxw(N1);  % finite difference points for 1st variable
    [y1,dy1]    = Maps.PhysSpace1(x1); %get physical space [0 L] 
    
    [x2,wInt2]  = FDxw(N2);  % finite difference points for 2nd variable
    [y2,dy2]    = Maps.PhysSpace2(x2); %get physical space [-half_wedge_angle half_wedge_angle]     
    
    % kron builds a large vector/matrix of the form ..    
  
    Pts     = struct('y1_kv',kron(y1,ones(size(y2))),...
                     'y2_kv',kron(ones(size(y1)),y2),...
                     'x1',x1,'x2',x2);    

    Diff     = GetDifferentiationMatrix(); 
    
    interpPts1=Maps.PhysSpace1( Pts.x1 );
    interpPts2=Maps.PhysSpace1( Pts.x2 );
    
    Interp   = FDFD_Interpolation(interpPts1,interpPts2,Pts,Maps);     
    
    Int      = GetIntegrationVector(); 
    %Ind      = GetIndicesBox(x1,x2);        
    Ind= [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Diff  = GetDifferentiationMatrix()      
    
    Diff1 = barychebdiff(x1,3);
    Diff2 = barychebdiff(x2,3);    
    
    Sel = {'Dy1'; 'DDy1'; 'DDDy1' ;'Dy2' ; 'DDy2' ; 'DDDy2'; ...
           'Dy1Dy2' ; 'DDy1Dy2' ; 'Dy1DDy2'; 'Lap' ; 'grad'}; 
    Diff = PhysicalDerivatives(Maps,Pts,Sel,Diff1,Diff2,3);    
      
    %First Variable:
    %[Dx1,DDx1,DDDx1]= barychebdiff(x1);        
    %[Dx2,DDx2,DDDx2]= barychebdiff(x2);            
    
    % NOTE THIS ONLY GIVES THE LAPLACIAN AND GRADIENT
    
    %Diff = PhysicalDerivatives(Maps,Pts,Sel,Dx1,Dx2,DDx1,DDx2,DDDx1,DDDx2);           
                  
end
 
function M  = GetIntegrationVector()
    M = kron(dy1'.*wInt1,dy2'.*wInt2);    
end



end
