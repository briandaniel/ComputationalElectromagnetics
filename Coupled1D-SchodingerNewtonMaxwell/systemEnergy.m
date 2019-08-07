function [Etot] = systemEnergy(param, N, Npsi, Nion, zDomHalf, dz, dt, maxStep, PSI, spin, Zion,  b1, b3, b4, gamma, e, hbar, me, epsFactor, fftCoul)
%SYSTEMENERGY Summary of this function goes here
%   Detailed explanation goes here

% Convert dt --> -i*dt for ground state evolution
dt = -1i*dt;

% mid = (zDomHalf(end) - zDomHalf(1))/2 + zDomHalf(1);
% coul = lorzScale.*gamma./sqrt(gamma^2 + ( zDomHalf - mid ).^2);
% fftcoul = fft([coul'; zeros(N-1,1)]);
% 
c1 = b1;
c2 = b4;
c3 = b3;

zMax = zDomHalf(end);
zMin = zDomHalf(1);
count = linspace(-1,1,Nion);
mid = (zDomHalf(end) - zDomHalf(1))/2 + zDomHalf(1);
Rion = mid + param*(zMax - zMin)./3*count; 

% Generate (n+1) time step matrix
u1 = ( -dt*b1/(2*dz^2)*ones(1,N-1) )                ;
u2 = ( (1+dt*b1/dz^2)*ones(1,N)  ) ;
u3 = ( -dt*b1/(2*dz^2)*ones(1,N-1) )                ;

% Generate (n) time step matrix
v1 = ( dt*b1/(2*dz^2)*ones(1,N-1) )                     ;
v2 = ( (1-dt*b1/(dz^2))*ones(1,N)  ) ;
v3 = ( dt*b1/(2*dz^2)*ones(1,N-1) )                     ;
B = diag(v1,-1) + diag(v2) + diag(v3, 1);

[Q R] = qr(PSI);
PSI = Q(:,1:Npsi);

for r = 1:maxStep
    Vx = zeros(1,N);
    for s = 1:Nion
    Vx = Vx - Zion(s)...
         *abs( 1./sqrt(gamma^2 + ( zDomHalf - Rion(s) ).^2) );      
    end
    Vx = Vx';
    % Calculate electron density for current time step
    rho = zeros(length(PSI),1);
    for s = 1:Npsi
       rho = rho + abs(PSI(:,s)).^2;
    end

    fftRho = fft( [rho; zeros(N-1,1)]);      
    Vcoul2  =  ifft(fftRho.* fftCoul); 
    Vcoul  = Vcoul2(floor(N/2+1):floor(3*N/2));
    
    %%%%%%%%%%% STRANG SPLITTING TO COMPPSIE NONLINEAR PART %%%%%%%%%%%%
    for s = 1:Npsi
        U = PSI(:,s);  
        U = exp(dt*(c3.*Vcoul)/2 + c2*dt*Vx/2).*U; 
        PSI(:,s) = U;
    end
    %%%%%%%%% END STRANG SPLITTING TO COMPPSIE NONLINEAR PART %%%%%%%%%%
    
    for s = 1:Npsi
        PSI(:,s) = B*PSI(:,s);
        a = u1; b = u2; c = u3;
        PSI(:,s) = tridisolve(a, b, c, PSI(:,s));
    end
    
    % Calculate electron density for current time step
    rho = zeros(length(PSI),1);
    for s = 1:Npsi
       rho = rho + abs(PSI(:,s)).^2;
    end

    fftRho = fft( [rho; zeros(N-1,1)]);      
    Vcoul2  =  ifft(fftRho.* fftCoul); 
    Vcoul  = Vcoul2(floor(N/2+1):floor(3*N/2));
    
    
    %%%%%%%%%%% STRANG SPLITTING TO COMPPSIE NONLINEAR PART %%%%%%%%%%%%
    for s = 1:Npsi
        U = PSI(:,s);  
        U = exp(dt*(c3.*Vcoul)/2 + c2*dt*Vx/2).*U; 
        PSI(:,s) = U;
    end
    %%%%%%%%% END STRANG SPLITTING TO COMPPSIE NONLINEAR PART %%%%%%%%%%

    
    % Continuous plotting to display solPSIion. Run time much faster if this
    % is removed.
    if mod(r,100) == 0;
        hold on
        plot(zDomHalf, real(PSI), zDomHalf, Vx/max(Vx), 'r--', zDomHalf, real(rho), 'b.' );
        plot(Rion, zeros(length(Rion),1), 'ro')
        text(-2, -.1, 'dt -> -i*dt ', 'fontsize', 18, 'fontname', 'cordia new')
        text(-2, -.14, strcat( 'time step: ', num2str(r)), 'fontsize', 18, 'fontname', 'cordia new')
        axis( [ zMin, zMax,-.3, .3] )
        getframe;
        hold off
        clf
    end
    
    
    for s = 1:Npsi
        PSI(:,s) = PSI(:,s)/norm(PSI(:,s));    
    end
        PSI = gramSchmidt(PSI);
%     % Normalize the solPSIion during evolPSIion to ground state.
%     PSIstring = zeros(N*Npsi/Zion(1), Zion(1));
%     for ss = 1: Zion(1)
%         PSItemp = [];
%         for s = 1:(Npsi/Zion(1))
%             PSItemp =  [PSItemp; PSI(:,Npsi/Zion(1)*(ss-1)+s)];
%         end
%         PSIstring(:,ss) = PSItemp;
%     end
% 
%     PSIstring = gramSchmidt(PSIstring);
%     size(PSIstring)
%     for ss = 1:Zion(1)
%         for s = 1:Npsi/Zion(1)
%             PSI(:, Npsi/Zion(1)*(ss-1)+s) =  PSIstring( (s-1)*N+1 : (s)*N, ss);
%         end
%     end
    
    
%     
%     PSI = gramSchmidt(PSI);
%     for s = 1:Npsi
%         PSI(:,s) = PSI(:,s)/sqrt(trapz(zDomHalf, abs(PSI(:,s)).^2));    
%     end
    
end

%% Calculate energy of system

% Second derivative of U using fourth order F.D. approximation.
%
%          -U(x+2h) + 16U(x+h) - 30U(x) + 16U(x-h) - U(x-2h)
% U_xx =  ---------------------------------------------------  + O(h^4)
%                              12h^2

% The corresponding matrix, D2 (Second derivative matrix)
e1 = -30*ones(1,N);
e2 = 16*ones(1,N-1);
e3 = -1*ones(1,N-2);
D2 = diag(e1) + diag(e2,-1) + diag(e2, 1) + diag(e3,-2) + diag(e3,2);
D2(end,1) = 16;
D2(end-1,1) = -1;
D2(end,2) = -1;
D2(1,end) = 16;
D2(1,end-1) = -1;
D2(2,end) = -1;
% Divide each element by constant 12*dx^2
D2 = D2./(12*dz^2);


% Electron-Electron repulsion potential
rho = zeros(length(PSI),1);
for s = 1:Npsi 
    rho = rho + (abs(PSI(:,s)).^2);
end
fftRho = fft( [rho; zeros(N-1,1)]);      
Vcoul2 =  ifft(fftRho.* fftCoul); 
Vcoul  = Vcoul2(floor(N/2+1):floor(3*N/2));    
Ue =.5*e^2*epsFactor*trapz( rho .* Vcoul ) ;  
if b3 == 0 
    Ue = 0;
end

% Electron kinetic energy
Te = 0;
for s = 1:Npsi 
    integrand = conj(PSI(:,s)).* ( D2* PSI(:,s) );
    Te = Te - spin*( hbar^2/(2*me) )*trapz(integrand);
end

% Electron-ion bonding potential
Ubond = 0;
for s = 1:Nion;
    integrand = rho'.*1./sqrt(gamma^2 + ( zDomHalf - Rion(s) ).^2);  
    Ubond = Ubond - spin*Zion(s)*e^2*epsFactor*trapz(integrand);
end


% Ion/Ion repulsion potential
Uion = 0;
for s = 1:Nion   
    for ss = 1:Nion
        if s ~= ss;
            const = .5*Zion(s)*Zion(ss)*e^2.*epsFactor;
            Enm   =  const./sqrt( gamma^2 + ( Rion(s) - Rion(ss) )^2 );
            Uion = Uion + Enm;
        end
    end
end

% dis = (Rion(2)-Rion(1))*1e9;
Etot = Te + Ue + Ubond + Uion;



end

