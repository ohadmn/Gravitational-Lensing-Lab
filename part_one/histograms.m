a_0 = 1;
b_0 = 2;
c_0 = 3;


a_val = [];
b_val = [];
c_val = [];
xi_val = [];

for i = 1:10000
    %[1,2,3,4,5]
    x = linspace(1,5,100);
    y_org = (a_0.*x.^2+b_0.*x+c_0);
    %y = (a.*x.^2+b.*x+c);
    y = y_org+0.001*randn(1,100);
    sigma = 0.1 * ones(1,100);
    %sigma = [0.1,0.1,0.1,0.1,0.1];
    [a,b,c] = parab(x,y,sigma);
    a_val(i)=a;
    b_val(i)=b;
    c_val(i)=c;
    xi_val(i) = sum((y-a.*x.^2-b.*x-c).^2./sigma.^2);

end

subplot(2, 2, 1);
histogram(a_val)
xlabel('it')
ylabel('a')
title('Histogram of a')
 
subplot(2, 2, 2);
histogram(b_val)
xlabel('it')
ylabel('b')
title('Histogram of b')
 
subplot(2, 2, 3);
histogram(c_val)
xlabel('it')
ylabel('c')
title('Histogram of c')

subplot(2, 2, 4);
histogram(xi_val)
xlabel('it')
ylabel('xi')
title('Histogram of xi')
