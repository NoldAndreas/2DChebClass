function A = GetTriangleParameters(Y)
%function A = GetTriangleParameters(Y)
%Input: Y(i,:) = (y_1,y_2), coordinates of i-th corner of triangle
%Output: 
% parameters such that
% y = A * (1;x_1;x_2;x_1*x_2)

    B = 1/4*[2,1,1;
            -2,1,1;
             0,-1,1;
             0,-1,1];
         
    A = (B*Y)';
end