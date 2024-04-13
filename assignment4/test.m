egf = 0.2*(1+(1/sqrt(3)))/2;
disp(egf);

ax = 0; bx = 0.2; cx = 1/sqrt(3);
tgf = gauss(ax, bx, cx);
disp(tgf);
% Function to compute the Gauss quadrature.
function result = gauss(a, b, c)
    % This function calculates the average of the square of a and b,
    % and the product of c and d.

    % Calculate the square of a and b, then find their average
    result = (a*(1 - c)/2) + (b*(1 + c)/2);
end