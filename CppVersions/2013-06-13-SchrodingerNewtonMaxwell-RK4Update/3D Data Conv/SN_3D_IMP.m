clear all
close all
fclose('all');

%% IMPORT DATA
folder = 'C:\Users\Inchoate\Documents\Visual Studio 2012\Projects\SNM update 5-15 - Copy\SNM\';
% folder = 'C:\Users\Brian\Documents\Visual Studio 2012\Projects\SNMdata\24elc24IonwInputfield\Run 3\SNM\';
% folder = 'C:\Users\Brian\Documents\Visual Studio 2012\Projects\batchRun\batchFiles\batch_0\dataSet_4\';

intData     = fopen( strcat(folder, 'intData.bin'       ) );
doubleData  = fopen( strcat(folder, 'doubleData.bin'    ) ); 
data        = fopen( strcat(folder, 'data.bin'          ) );
vidData     = fopen( strcat(folder, 'vid.bin'           ) ); 
cxData     = fopen( strcat(folder, 'vidCX.bin'           ) ); 

indexData = fread(intData, 'int');
d2Data = fread(doubleData, 'double');
vid = fread(vidData, 'double');
vidCX = fread(cxData, 'double');

%------------------------------------------------------------------------%

N3 = indexData(1);
Nx = indexData(2);
Ny = indexData(3);
Nz = indexData(4);
frames = indexData(5);

data = fread(data, 'double');

x = linspace(0,10e-9,Nx);
dx = x(2) - x(1);


Uvec = data(1:N3);
Vvec = data(N3+1:end);

U = reshape(Uvec,Nz,Ny,Nx);
Vext = reshape(Vvec,Nz,Ny,Nx);

rho = U(:,:,25);
M(:,:) = Vext(25,:,:);

%------------------------------------------------------------------------%

figure('outerposition', [100 140 600 1000]);
subplot(2,1,1)
    surf(x, x, M);  

    rhoMax = max(max(max(rho)));
subplot(2,1,2)
   surf(x, x, rho/rhoMax);  


K = jet(4);
A = figure('outerposition', [10 100 1800 1050]);

vidMax = max(abs(vid(50*N3:end)));
uMax = max(abs(vidCX(50*N3:end)));

%------------------------------------------------------------------------%
        
for k = 0:3:2*frames-1
 
    clf
    %---------------------------------------------------------------------%
        
    rhoVec = vid(k*N3+1:(k+1)*N3);
    rho = reshape(rhoVec,Nz,Ny,Nx);
    
    uVec = vidCX(k*N3*2+1:(k+1)*N3*2);
    U = uVec(1:2:end-1) + 1i*uVec(2:2:end);
    U = reshape(U,Nz,Ny,Nx);
    
    %---------------------------------------------------------------------%
        
    subplot(2,3,1)
        plot(x,  rho(:,Ny/2,Nx/2));
        text(x(end/2), .6 , num2str(k) )
        axis([x(1) x(end) -vidMax vidMax ])
  
    text(0,0,num2str(k));

    subplot(2,3,2)
        surf(x, x, rho(:,:,Nx/2), 'edgealpha', .3);  
        axis([x(1) x(end) x(1) x(end) -vidMax  vidMax ])
            

    subplot(2,3,3)
        p1 = patch( isosurface(x,x,x,rho,.01));
        p2 = patch( isosurface(x,x,x,rho,-.01));
        
        set(p1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'facealpha', 1);
        set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'facealpha', 1);
        isonormals(x,x,x, rho, p1)
        isonormals(x,x,x, rho, p2)

        daspect([1,1,1])
        view(3); 
        axis([x(1) x(end) x(1) x(end) x(1) x(end)])
        camlight 
        lighting gouraud
    %---------------------------------------------------------------------%
        
    subplot(2,3,4)
        plot(x, real(U(:,Ny/2,Nx/2)));
        text(x(end/2), .6 , num2str(k) )
        axis([x(1) x(end) -uMax  uMax ])
  
    text(0,0,num2str(k));

    subplot(2,3,5)
        surf(x, x,real(U(:,:,Nx/2)), 'edgealpha', .3);  
        axis([x(1) x(end) x(1) x(end) -uMax  uMax ])
            

    subplot(2,3,6)
        p1 = patch( isosurface(x,x,x,real(U),uMax/4));
        p2 = patch( isosurface(x,x,x,real(U),-uMax/4));
        
        set(p1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'facealpha', 1);
        set(p2, 'FaceColor', 'red', 'EdgeColor', 'none', 'facealpha', 1);
        isonormals(x,x,x, rho, p1)
        isonormals(x,x,x, rho, p2)

        daspect([1,1,1])
        view(3); 
        axis([x(1) x(end) x(1) x(end) x(1) x(end)])
        camlight 
        lighting gouraud
   %---------------------------------------------------------------------%
        
    getframe;
    
end







