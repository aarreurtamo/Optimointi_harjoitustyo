%% BM20A8200 Harjoitustyö 4
clearvars
close all
clc

% Minimoitava funktio
%f = @(x) -z(x(1),x(2));
f = @(x) 1e-2*x(1).^2+1e-2*x(2).^2+7*1e-3*x(1).*x(2)...
    -144*x(1)-114*x(2)+4*1e5;
% Alkuarvaus
x0 = [0,0]';

% fmincon
A = [1,0;...
    0,1;...
    1,1];
b = [5*10^3,8*10^3,10^4];
lb = [0,0];
ub = [];
x = fmincon(f,x0,A,b,[],[],lb,ub,[]);

disp('Minimipiste fmincon:')
disp(x)


%% Sakkomenetelmä 
clearvars
close all
clc

% Minimoitava funktio.
f = @(x) 1e-2*x(1).^2+1e-2*x(2).^2+7*1e-3*x(1).*x(2)...
    -144*x(1)-114*x(2)+4*1e5;

% Alkuarvaus
x0 = [0,0]';

% Toleranssi
tol = 1e-6;

% Rajoittamaton 
g1 = @(x) -[x(1);x(2)];
[x_min1,ite1,sakkokierrokset1] = ngs(f,g1,x0,tol);

% Rajoitettu
g2 = @(x) [x(1)-5*1e+3;...
          x(2)-8*1e+3;...
          x(1)+x(2)-10^4;...
          -x(1);...
          -x(2)];
[x_min2,ite2,sakkokierrokset2] = ngs(f,g2,x0,tol);



disp('Minimipiste (rajaoittamaton kapasiteetti): ')
disp(x_min1)
disp('Minimipiste (rajoitettu kapasiteetti): ')
disp(x_min2)


%% Defiinttisyyden tarkastelu

Hf = [2*1e-2,7*1e-3;...
        7*1e-3,2*1e-2];

omin = eig(Hf)


%% Piirtokomennot

% Maksimoitava voittofuktio
z = @(x1,x2) -1e-2*x1.^2-1e-2*x2.^2-7*1e-3*x1.*x2+144*x1+114*x2-4*1e+5;

Raja = 1e+4;
x1x1 = linspace(0,Raja,50);
x2x2 = linspace(0,Raja,50);
[X1,X2] = meshgrid(x1x1,x2x2);

Z = z(X1,X2);

figure;
mesh(X1,X2,Z)
hold on
plot3(x_min1(1),x_min1(2),z(x_min1(1),x_min1(2)),'k*',MarkerSize=10)

xlabel('x1 (Laadun 24 määrä)')
ylabel('x2 (Laadun 27 määrä)')
zlabel('z (voitto)')
title('Tilanne ilman rajallista valmistuskapasiteettia')
legend('Voittofunktio z','Maksimipiste')


% Piirtokomennot
x1x1 = linspace(0,5*1e+3,50);
x2x2 = linspace(0,8*1e+3,50);
[X1,X2] = meshgrid(x1x1,x2x2);


z1 = @(x1,x2) z(x1,x2).*(x1+x2-10^4 <= 0);
Z = z1(X1,X2);

figure;
mesh(X1,X2,Z)
hold on
plot3(x_min2(1),x_min2(2),z(x_min2(1),x_min2(2)),'k*',MarkerSize=10)

xlabel('x1 (Laadun 24 määrä)')
ylabel('x2 (Laadun 27 määrä)')
zlabel('z (voitto)')
title('Tilanne rajallisella valmistuskapasiteetilla')
legend('Voittofunktio z','Maksimipiste')

