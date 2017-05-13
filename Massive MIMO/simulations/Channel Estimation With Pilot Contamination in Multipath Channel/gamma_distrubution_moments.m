clear all;close all;clc

M = 10;
theta_ik = 1.3;

N = 100000;
z = zeros(1,N);
for i=1:1:N
    y = 0;
    for m=1:1:M
        x = (1/sqrt(2))*complex(randn(1,1),randn(1,1));
        y = y + theta_ik*abs(x).^2;
    end
    
    z(i) = y;
end

fprintf(1,'Estimated mean:¨%d\n',mean(z));
fprintf(1,'Actual mean:¨%d\n\n\n\n\n',M*theta_ik);

%
fprintf(1,'Estimated variance:¨%d\n',var(z));
fprintf(1,'Actual variance:¨%d\n',M*(theta_ik.^2));