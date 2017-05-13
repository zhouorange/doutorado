clear all;close all;clc

M = 10;
P = 4;
epsilon_ik = 1.3;

N = 1000000;
z = zeros(1,N);
for i=1:1:N
    q = 0;
    for m=1:1:M
        y0 = sqrt(epsilon_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        y1 = sqrt(epsilon_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        
        q = q + y0*y1;        
    end
    
    z(i) = q;
end

fprintf(1,'Estimated mean:¨%d\n',mean((z)));
fprintf(1,'Estimated var:¨%d\n',var((z)));

figure;hist(real(z),100)
figure;hist(imag(z),100)
