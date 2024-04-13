% Assuming densities and sensitivities are loaded or defined here
nelx = 1; nely = 20; volfrac = 0.3;

% Call the modified OC function
xnew = OC(nelx, nely, density, volfrac, sensitivity);

% Prepare the output matrix
outputMatrix = [density(:), sensitivity(:), xnew(:), xnew(:) - density(:)];

% Display the output matrix
disp('Densities, Sensitivities, NewDensities, Differences:');
for i = 1:size(outputMatrix, 1)
    fprintf('%10.4f %10.4f %10.4f %10.4f\n', outputMatrix(i, 1), outputMatrix(i, 2), outputMatrix(i, 3), outputMatrix(i, 4));
end

function [xnew]=OC(nelx, nely, x, volfrac, dc)
    l1 = 0; l2 = 100000; move = 0.25;
    while (l2-l1 > 1e-4)
        lmid = 0.5*(l2+l1);
        xnew = max(0.001, max(x-move, min(1., min(x+move, x.*((-dc./lmid).^0.6)))));
        if sum(xnew(:)) - volfrac*nelx*nely > 0
            l1 = lmid;
        else
            l2 = lmid;
        end
    end
end
