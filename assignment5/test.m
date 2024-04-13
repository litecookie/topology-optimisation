f = @(x) (exp(x / 200) * (sin((pi * x) / 500) + 2)) - 25;
x_value = 514.375738;
result = f(x_value);

disp(['The result of substituting x = ', num2str(x_value), ' into the function is: ', num2str(result)]);
