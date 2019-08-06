clear all
close all


intData = fopen('C:\Users\Brian\Documents\Visual Studio 2012\Projects\PoissonSolver\poissonMG\intData.bin');
datas = fopen('C:\Users\Brian\Documents\Visual Studio 2012\Projects\PoissonSolver\poissonMG\data.bin');

indexData = fread(intData, 'int');
N3 = indexData(1);
Nx = indexData(2);
Ny = indexData(3);
Nz = indexData(4);

data = fread(datas, 'double');

Uvec = data(1:N3);
fvec = data(N3+1:2*N3);
UexactVec = data(2*N3+1:3*N3);
guessVec = data(3*N3+1:4*N3);
dUvec = data(4*N3+1:5*N3);

U = reshape(Uvec,Nz,Ny,Nx);
f = reshape(fvec,Nz,Ny,Nx);
Uexact = reshape(UexactVec, Nz, Ny, Nx);
guess = reshape( guessVec, Nz, Ny, Nx );
dU = reshape( dUvec, Nz, Ny, Nx);

xMax = 6.0;
yMax = 6.0;
zMax = 6.0;

x = linspace(0,xMax,Nx);
y = linspace(0,yMax,Ny);
z = linspace(0,zMax,Nz);

dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);
% 
% [X , Y] = meshgrid(x,y);
% 
% Z = ( sin(X*6) ) .* exp( - (X - 3).^2 -  (Y - 3).^2);
% 
% surf(X,Y, Z, 'edgealpha', .1 );
% 
% 
% return

%-------------------------------------------------------------------------%
M = zeros(Ny,Nx);
Mexc = zeros(Ny,Nx);
Mg = zeros(Ny,Nx);
Mf = zeros(Ny,Nx);
Md = zeros(Ny,Nx);

M(:,:) = U(round(Nz/2.0),:,:);
Mexc(:,:) = Uexact(round(Nz/2.0),:,:); 
Mg(:,:) = guess(round(Nz/2.0),:,:); 
Md(:,:) = dU(round(Nz/2.0),:,: ); 
Mf(:,:) = f(round(Nz/2.0),:,: );


figure('outerposition', [100 100 1600 1000]);

subplot(2,3,1)
    surf(y, x, Mf', 'edgealpha', .1);   
title ( '\Delta U = f; f' )

subplot(2,3,2)
    surf(y, x, M', 'edgealpha', .1);  
title ( '\Delta U = f; Uapprox' )
axis([0 yMax 0 xMax -1.1 1.1])

subplot(2,3,3)
    surf(y, x, Mexc', 'edgealpha', .1);
title ( '\Delta U = f; Uexact' )
axis([0 yMax 0 xMax -1.1 1.1])

subplot(2,3,4)
    surf(y, x, Mg', 'edgealpha', .1);    
title ( '\Delta U = f; guess for U' )

subplot(2,3,5)
    surf(y, x, M'-Mexc', 'edgealpha', .1);    
title ( 'Error (U - Uexact)' )


%-------------------------------------------------------------------------%
M = zeros(Nz,Ny);
Mexc = zeros(Nz,Ny);
Mg = zeros(Nz,Ny);
Mf = zeros(Nz,Ny);
Md = zeros(Nz,Ny);

M(:,:) = U(:,:,round(Nx/2.0));
Mexc(:,:) = Uexact(:,:,round(Nx/2.0)); 
Mg(:,:) = guess(:,:,round(Nx/2.0)); 
Md(:,:) = dU(:,:,round(Nx/2.0)); 
Mf(:,:) = f(:,:,round(Nx/2.0));


figure('outerposition', [100 100 1600 1000]);

subplot(2,3,1)
    surf(z, y, Mf', 'edgealpha', .1);   
title ( '\Delta U = f; f' )

subplot(2,3,2)
    surf(z, y, M', 'edgealpha', .1);  
title ( '\Delta U = f; Uapprox' )
axis( [z(1) z(end) y(1) y(end) 0 1.1])

subplot(2,3,3)
    surf(z, y, Mexc', 'edgealpha', .1);
title ( '\Delta U = f; Uexact' )
axis( [z(1) z(end) y(1) y(end) 0 1.1])

subplot(2,3,4)
    surf(z, y, Mg', 'edgealpha', .1);    
title ( '\Delta U = f; guess for U' )

subplot(2,3,5)
    surf(z, y, M'-Mexc', 'edgealpha', .1);    
title ( 'Error (U - Uexact)' )
%-------------------------------------------------------------------------%






