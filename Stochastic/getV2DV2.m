function [V2,DV2]=getV2DV2(x,optsPhys)

switch optsPhys.type
    
    case 'stoc'
        
        x=x(:);
        
        % input is the full vector of coordinates so get the separations
        dim=optsPhys.dim;
        nParticles=length(x)/dim;
        Rij=getRij(x,x,dim);
      
    otherwise
        
        % input is already interparticle distances
        Rij=x;
        
end

% get name of function
V2DV2=str2func(optsPhys.V2DV2);

% and evaluate it.  Note in particular that the derivative already contains
% the 1/r from from grad|x| = x/|x|.
[V2,V2prime_r]=V2DV2(Rij,optsPhys);

switch optsPhys.type
    
    case 'stoc'

        %------------------------------------------------------------------
        % Calculate V2                                                 
        %------------------------------------------------------------------

        % factor of 1/2 to avoid double counting and remove diagonal
        % elements as there should be no self-interaction in the stochastic
        % calculations (although this doesn't actually change anything as
        % it's equivalent to shifting the external potential by a constant)
        V2 = V2-diag(diag(V2));       
        V2 = 1/2*sum(sum(V2,2),1);  

        %------------------------------------------------------------------
        % Calculate DV2                                                        
        %------------------------------------------------------------------

        % each column will contain the derivative in that coordinate
        DV2=zeros(nParticles,dim);

        % rows of x now contain the dth coordinate of each particle
        x=reshape(x,dim,nParticles);

        % remove diagonal elements as no self-interaction
        V2prime_r = V2prime_r -diag(diag(V2prime_r));
        
        for d=1:dim
            % calculate the derivative for each coordinate
            % the i,j element of bsxfun(@minus,(x(d,:))',x(d,:)) is 
            % x_(i,d)-x_(j,d), the difference in the dth coordinate of
            % particles i and j.
            % We then elementwise multiply by the derivative of V2 (including
            % the 1/|x| from grad|x| = x/|x|) and then sum over the second
            % (column) coordinate which corresponds to summing over j.  Hence
            % the ith row and dth column of DV2 is DV2 for coordinate d of
            % particle i

            DV2(:,d)=sum(V2prime_r.*bsxfun(@minus,(x(d,:))',x(d,:)),2);

        end

        % take the transpose so that there are dim rows containing the
        % gradients for each particle.  Then form into a single vector.
        DV2=DV2';
        DV2=DV2(:);
    
    otherwise  % don't need this for ddft or plotting
        
        DV2=[];
end