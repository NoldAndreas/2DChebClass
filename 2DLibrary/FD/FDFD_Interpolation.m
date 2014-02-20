function InterpValues = FDFD_Interpolation(interp1,interp2,xPts,Maps)
    
    % Doesn't do any interpolation: if the input points are not the same as
    % the physical points then throw an error.
    
    
    
    
    
    if( isequal(interp1 , Maps.PhysSpace1( xPts.x1 ) ) && ...
            isequal(interp2 , Maps.PhysSpace2( xPts.x2 ) ) )
        
        Nplot1 = length(interp1);
        Nplot2 = length(interp2);

        % can use identity as the vectors we apply it to just store the
        % values and it doesn't matter whether we plot in the computational
        % or physical space.
        Interp1 = speye(Nplot1); 
        Interp2 = speye(Nplot2);
    
        InterPol = sparse( kron(Interp1,Interp2) );
        pts1 = sparse( kron(interp1,ones(size(interp2))) );
        pts2 = sparse( kron(ones(size(interp1)),interp2) );
        
        % kron form of the two interpolations                                                            
        InterpValues = struct('InterPol', InterPol,...
                            'pts1', pts1,...
                            'pts2', pts2,...
                            'Nplot1',Nplot1,...
                            'Nplot2',Nplot2);
                        
    else
        error('Input vectors should agree with physical points');
    end
    
end