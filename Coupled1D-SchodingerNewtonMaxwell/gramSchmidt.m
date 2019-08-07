function [ Q ] = gramSchmidt(UT )
%GRAMSCHMIDT Computes orthonormal columns according to the Gram-Schmidt
%orthonormalization procedure.
    
    sU = size(UT);
    Q = zeros(sU);    
    for k = 1:sU(2)
        v = UT(:,k);
        for n = 1:k-1;
            R = Q(:,n)'*UT(:,k);
            v = v - R*Q(:,n);
        end
 
        Q(:,k) = v/norm(v);   

    end
   
end

