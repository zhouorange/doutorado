function [ft] = calculateFt(t, M)

ft = (gamma(2*M)/(gamma(M)^2))*((t*(1-t))^(M-1));
