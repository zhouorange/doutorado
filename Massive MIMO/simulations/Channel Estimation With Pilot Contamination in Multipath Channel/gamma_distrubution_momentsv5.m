clear all;close all;clc

M = 10;
P = 4;
theta_ik = 1.3;

N = 10000000;
z = zeros(1,N);
for i=1:1:N
    q1 = 0;
    q2 = 0;
    q3 = 0;
    q4 = 0;
    q5 = 0;
    for m=1:1:M
        y1 = sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        q1 = q1 + abs(y1).^2;
        
        y2 = sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        q2 = q2 + abs(y2).^2;
        
        y3 = sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        q3 = q3 + abs(y3).^2;
        
        y4 = sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1));
        q4 = q4 + abs(y4).^2;
        
        q5 = q5 + (y1 * sqrt(theta_ik)*(1/sqrt(2))*complex(randn(1,1),randn(1,1)));
        
    end
    
    z(i) = q5 ./ ((q1 + q2 + q3 + q4).^2);
end

fprintf(1,'Estimated mean:¨%d\n',mean(z));
%fprintf(1,'Actual mean:¨%d\n\n\n\n\n', 1./(theta_ik*P*(P*M-1)) );
