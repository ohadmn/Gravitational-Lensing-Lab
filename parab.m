function [a,b,c] = parab(x,y,sigma)
    a_11 = sum(x.^4./sigma.^2);
    a_12 = sum(x.^3./sigma.^2);
    a_13 = sum(x.^2./sigma.^2);
    a_21 = sum(x.^3./sigma.^2);
    a_22 = sum(x.^2./sigma.^2);
    a_23 = sum(x.^1./sigma.^2);
    a_31 = sum(x.^2./sigma.^2);
    a_32 = sum(x.^1./sigma.^2);
    a_33 = sum(x.^0./sigma.^2);

    mat = [a_11 a_12 a_13;
           a_21 a_22 a_23;
           a_31 a_32 a_33];

    b_1 = sum(y.*x.^2./sigma.^2);
    b_2 = sum(y.*x.^1./sigma.^2);
    b_3 = sum(y.*x.^0./sigma.^2);

    RH_vec = [b_1, b_2, b_3];
    inv_mat = inv(mat);

    %parab = inv_mat * RH_vec;
    parab = RH_vec/mat;
    
    a = parab(1);
    b = parab(2);
    c = parab(3);

end