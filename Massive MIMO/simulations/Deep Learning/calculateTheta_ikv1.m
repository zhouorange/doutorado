function [theta] = calculateTheta_ikv1(beta_iik, epsilon11, M)

kik = calculateKik(beta_iik, epsilon11);

factor1 = (gamma(2.*M)./(gamma(M).^2));

factor2 = exp(gammaln(2.*M)-2.*gammaln(M));

fun = @(t,w) ((((kik.^2).*(1-t)+kik.*w.*sqrt(t.*(1-t)))./((kik.^2).*(1-t)+2.*kik.*w.*sqrt(t.*(1-t))+t)).*((gamma(2.*M)./(gamma(M).^2)).*((t.*(1-t)).^(M-1))).*((M./pi).*beta(0.5,M).*((1-(w.^2)).^(M-0.5))));

theta = integral2(fun,0,1,-1,1,'Method','auto','AbsTol',1e-12,'RelTol',1e-10);
%theta = integral2(fun,0,1,-1,1);

