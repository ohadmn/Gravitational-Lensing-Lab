x_val = linspace(-62,62,104);
y_val = [];
sigma_0 = linspace(1,1,104);
disp(sigma_0);
for i = 1:104
    y_val(i) = (-5*x_val(i)^2-4*x_val(i)+50)+10*randn(1);

end

[a,b,c] = parab(x_val,y_val,sigma_0);

plot(y_val)

