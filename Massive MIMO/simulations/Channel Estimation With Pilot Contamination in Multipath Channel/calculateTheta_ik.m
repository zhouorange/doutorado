function [theta] = calculateTheta_ik(beta_iik, epsilon11, M)

t_step = 0.01;
w_step = 0.01;

kik = calculateKik(beta_iik, epsilon11);

theta = 0;

for t=0:t_step:1
    
    ft = calculateFt(t, M);
    
    for w=-1:w_step:1
       
        fw = calculateFw(w, M);
            
        aux = (((kik^2)*(1-t)+kik*w*sqrt(t*(1-t))) / ((kik^2)*(1-t)+2*kik*w*sqrt(t*(1-t))+t)) * ft * fw;
        
        theta = theta + aux;
        
    end
end

%numOfIter = ((1/t_step) + 1)*((1/w_step)*2 + 1);

%theta = theta/numOfIter;
