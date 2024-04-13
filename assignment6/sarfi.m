clc 
clear all

load("problem_set_06_data.mat")
 
%Given Data
E= 70000; % in N/mm2
nu= 0.3;
volfrac=0.3;
penal= 3;
ystrength= 3
85; % in N/mm2
gam = 0.2;
[nelx, nely] = size(rho);
% meshing 
x = [1:nelx]-0.5; % changing co ordinates to the centre of the element 
y = -(x);

[X,Y] = meshgrid(x,y);
 stressvm= zeros(nelx,nely); 
        A = 1;
for elx = 1:nelx
    for ely= 1:nely
       
        D = (rho(elx,ely)^penal)*(E/(1-nu^2))[1 nu 0;nu 1 0;0 0 (1-nu)/2];

        d = [dx(ely+1,elx) dy(ely+1,elx); dx(ely+1,elx+1) dy(ely+1,elx+1); dx(ely,elx+1) dy(ely,elx+1); dx(ely,elx) dy(ely,elx)];

        strain = 0.5*[d(2,1)+d(3,1)-d(4,1)-d(1,1); 
                      d(3,2)+d(4,2)-d(1,2)-d(2,2); 
                      d(2,2)+d(3,1)+d(3,2)+d(4,1)-d(1,1)-d(1,2)-d(2,1)-d(4,2)];

        stress = D*strain; 
        sigmaxx = stress(1); 
        sigmayy = stress(2); 
        sigmaxy = stress(3);
        stressvm(ely,elx) = sqrt((sigmaxx.^2 + sigmayy.^2) + 3*sigmaxy.^2 - sigmaxx.*sigmayy);
    end
end

adjVm = stressvm./(rho.^penal);
relVm = zeros(ely,elx);
for elx =1:nelx
    for ely = 1:nely
        relVm(ely,elx) = stressvm(ely,elx)/((rho(ely,elx)^penal)*(1-gam+(gam/rho(ely,elx))));
    end
end

%contourf(X,Y,stress_vm,100,'edgecolor','none')
figure(1)
contourf(X,Y,adjVm,100,'edgecolor','none')
title('Adjusted Von-Misses Stress')
colorbar()

figure(2)
contourf(X,Y,relVm,100,'edgecolor','none')
title('Relaxed Von-Mises Stress')
colorbar()