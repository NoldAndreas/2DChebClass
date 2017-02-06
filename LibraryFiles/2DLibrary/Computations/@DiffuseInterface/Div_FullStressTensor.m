 function [A,b] = Div_FullStressTensor(this,phi)

    Diff      = this.IDC.Diff;    

    [A11,b11] = FullStressTensorIJ(this,phi,1,1); 
    [A12,b12] = FullStressTensorIJ(this,phi,1,2); 
    [A21,b21] = FullStressTensorIJ(this,phi,2,1); 
    [A22,b22] = FullStressTensorIJ(this,phi,2,2); 

    A1        = Diff.Dy1 * A11  + Diff.Dy2*A21;
    A2        = Diff.Dy1 * A12  + Diff.Dy2*A22;

    b1        = Diff.Dy1 * b11  + Diff.Dy2*b21;
    b2        = Diff.Dy1 * b12  + Diff.Dy2*b22;

    A   = [A1;A2];
    b   = [b1;b2];      

end