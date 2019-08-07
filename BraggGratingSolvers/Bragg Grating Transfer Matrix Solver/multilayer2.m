na = 1; nb = 1.52; nH = 2.32; nL = 1.00;        % refractive indices
LH = 0.25; LL = 0.25;                           % optical thicknesses in units of λ0


la0 = 500;                                      % λ0 in units of nm
rho = (nH-nL)/(nH+nL);                          % reflection coefficient ρ
la2 = pi*(LL+LH)*1/acos(rho) * la0;             % Right bandedge
la1 = pi*(LL+LH)*1/acos(-rho) * la0;            % left bandedge
Dla = la2-la1;                                  % bandwidth

N = 8;                                          
n = [na, nH, repmat([nL,nH], 1, N), nb];        % indices for the layers A|H(LH)N |G 
L = [LH, repmat([LL,LH], 1, N)];                % Lengths of layers
la = linspace(300,800,501);                     % plotting range is 300 ≤ λ ≤ 800 nm
Gla = 100*abs(multidiel(n,L,la/la0)).^2; 


figure; plot(la,Gla); 
f = linspace(0,6,1201); 
Gf = 100*abs(multidiel(n,L,1./f)).^2; 
figure; plot(f,Gf); 
