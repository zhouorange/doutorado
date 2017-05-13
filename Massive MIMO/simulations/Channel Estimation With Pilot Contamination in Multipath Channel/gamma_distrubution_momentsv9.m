clear all;close all;clc

M = 30;
P = 4;
theta_ik = 1.3;

N = 1000000;
z = zeros(1,N);
for i=1:1:N
    q1 = 0;
    q2 = 0;
    q3 = 0;
    q4 = 0;
    for m=1:1:M
        y1 = sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        q1 = q1 + abs(y1).^2;
        
        y2 = sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        q2 = q2 + abs(y2).^2;
        
        y3 = sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        q3 = q3 + abs(y3).^2;
        
        y4 = sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        q4 = q4 + abs(y4).^2;
    end
    
    z(i) = q1 + q2 + q3 + q4;
end

fprintf(1,'Estimated mean:¨%d\n',mean(z));
fprintf(1,'Expected mean: %d\n\n\n',M*P*theta_ik);

fprintf(1,'Estimated variance:¨%d\n',var(z));
fprintf(1,'Expected variance:¨%d\n',M*P*theta_ik*theta_ik);

% figure;hist(real(z),1000)
% figure;hist(imag(z),1000)
