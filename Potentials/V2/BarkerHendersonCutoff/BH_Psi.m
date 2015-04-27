function z = BH_Psi(y,V2)
    
    rc  = V2.r_cutoff;        

    markG1 = (y > 1) & (y < rc);
    markL1 = (y <= 1);

    yG1 = y(markG1);
    yL1 = y(markL1);

    z  = zeros(size(y));
    z(markG1) = -0.2e1 / 0.45e2 * pi * (45 * yG1.^ 10 * rc^ 6 - 60 * yG1.^ 9 * rc ^ 7 + 15 * yG1.^ 6 * rc ^ 10 - 18 * yG1.^ 10 + 20 * yG1.^ 9 * rc - 2 * rc ^ 10) / (rc ^ 10)./(yG1.^ 9);
    z(markL1) = 0.2e1 * pi * (0.1e1 / rc ^ 4 - 0.2e1 / 0.5e1 / rc ^ 10) * (1 - yL1) - 0.6e1 / 0.5e1 * pi * (1 - yL1) - 0.2e1 / 0.45e2 * pi * (0.13e2 * rc ^ 10 - 0.60e2 * rc ^ 7 + 0.45e2 * rc ^ 6 + 0.20e2 * rc - 0.18e2) / rc ^ 10;


    %Rescale
    alpha     = 1/2*( BH_0D_I(rc) - BH_0D_I(1) );
    c         = V2.epsilon*(-16/9*pi)/alpha;    
    z         = c*z;    
    alpha     = c*alpha;
end