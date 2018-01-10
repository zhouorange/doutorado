function [kik] = calculateKik(beta_iik, epsilon11)

kik = sqrt(beta_iik/(epsilon11 - beta_iik));
