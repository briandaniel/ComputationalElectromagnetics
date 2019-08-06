function [ dUT ] = nlsSimple( t, UT, k, Vx, c1, c2) %#ok<INUSL>
% Solves the NLS using a pseudospectral fft method
% input requires variables and parameters 
%     ( t, UT, k, Vx, fftCoul, c1, c2, c3)
% i.e. time dependent schrodinger equation with linear potential and 
% nonlinear square term...
%-------------------------------------------------------------------------%
%                                                                         %
%        dPsi_i      d^2Psi_i                                             %
%        ------ = c1*-------- + c2*V(x)*Psi_i  + c3*F(rho(x,t))           %
%          dt          dx^2                                               %
%                                                                         %
%-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -%
%      where rho(x) = sum |Psi_i(x)|^2                                    %
%                                                                         %
%      F(rho) = Int [ rho(x)*u(x-x') ] dx'                                %
%                i.e. F is the convolution of rho(x) and u(x)             %
%                                                                         %
%-------------------------------------------------------------------------%

    pkg = size(UT);
    rho = zeros(pkg(1),1);
    dUT = zeros(pkg);
    for s = 1:pkg(2) 
        rho = rho + (abs(UT(:,s)).^2);
    end

    for s = 1:pkg(2)
        uxx = (ifft( (k*1i).^2 .* fft(UT(:,s)) ) ); 
        dUT(:,s) = c1.*uxx + ( c2*Vx ).*UT(:,s);
    end

end

