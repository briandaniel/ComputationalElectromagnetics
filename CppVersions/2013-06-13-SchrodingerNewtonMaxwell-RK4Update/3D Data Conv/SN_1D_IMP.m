clear all
close all
fclose('all');

%% IMPORT DATA
folder = 'C:\Users\Inchoate\Documents\Visual Studio 2012\Projects\SNM update 5-15 - Copy\SNM\';
% folder = 'C:\Users\Brian\Documents\Visual Studio 2012\Projects\SNMdata\24elc24IonwInputfield\Run 3\SNM\';
% folder = 'C:\Users\Brian\Documents\Visual Studio 2012\Projects\batchRun\batchFiles\batch_0\dataSet_4\';

data = fopen( strcat(folder, '1Ddata.bin' ) );
vidData = fopen( strcat(folder, '1Dvid.bin' ) );

data = fread(data, 'double');
vidData = fread(vidData, 'double');

N = data(1);
L = data(2);
dx = data(3);
x = 0:dx:L;


data = data(4:end);
U = data(1:N);
dU = data(N+1:end);

plot( x, U, x, dU)


frames = 200000;
figure;
for k = 0:1:frames;
   
    rho = vidData(k*N+1: (k+1)*(N));
    
    plot(x,rho);
     axis ([x(1), x(end), -1,1])
    text(0,0,num2str(k));
    getframe;
   
    
end
    
    
        
    
