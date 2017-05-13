function [tap_value, tap_pos, res]=mp(A,y)

[Np Ncp] = size(A);

Al = A;

for jj=1:1:Ncp
    norm = (sqrt(sum((abs(A(:,jj)).^2))));
    A(:,jj) = A(:,jj)./norm;
end

residual = y;

norm2_res = sqrt(sum(abs(residual).^2));

tap_value = complex(zeros(1,Ncp),zeros(1,Ncp));
tap_pos = zeros(1,Ncp);
res = complex(zeros(1,Ncp),zeros(1,Ncp));

for ii=1:1:Ncp
    
    scalarproducts = (A'*residual);
    
    [v,p] = max((abs(scalarproducts).^2));
    
    Theta_p = A(:,p);
    
    norm2 = sqrt(sum(abs(Theta_p).^2));
    
    %tap(ii) = scalarproducts(p)/norm;
    tap_value(ii) = Al(:,p)'*residual./(sum(Al(:,p)'*Al(:,p)));
    tap_pos(ii) = p;
    res(ii) = norm2_res;
    
    residual = residual - ((scalarproducts(p)*Theta_p)/norm2);
    
    A(:,p) = complex(zeros(Np,1),zeros(Np,1));
    
    norm2_res = sqrt(sum(abs(residual).^2));
    
end


