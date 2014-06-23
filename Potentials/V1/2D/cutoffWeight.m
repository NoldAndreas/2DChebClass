function [weight,Dweight] = cutoffWeight(r,Rin,Rout)

    R = (r-Rin)/(Rout-Rin);

    weight = 1 - exp(-(R./(1-R)).^2);
    weight(r<Rin) = 0;
    weight(r>Rout) = 1;

    Dweight = exp(-(R./(1-R)).^2).*(2*R.*(1-R).^(-3))/(Rout-Rin);
    Dweight(r<Rin | r>Rout) = 0;

end