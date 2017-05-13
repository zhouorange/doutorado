function [fw] = calculateFw(w, M)

fw = (M/pi)*beta(1/2,M)*((1-(w^2))^(M-0.5));