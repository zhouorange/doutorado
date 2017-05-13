function y = A_Fourier(z,pil_pos,D,J,mat,mode)

%%%     Copyright: (c) Georg Tauböck, 2006-2010
%%%
%%%     Generates CS Measurement Matrix

pilot_positions=pil_pos(:,mat);
s2=size(z,2);

if mode==1
    
    for kk=1:s2
        a(:,kk)=vec((repmat((exp((-1j)*pi*((0:(J-1))'))),1,D).*(ifft(((sqrt(J)/sqrt(D))*fft((vec2mat(z(:,kk),D)).')).'))).');
    end
    
    y=a(pilot_positions,:)*sqrt(J*D/length(pilot_positions));
    
else if mode==2
        
        a=z(1:(length(pilot_positions)),:);
        b=zeros(D*J,s2);
        b(pilot_positions,:)=a;
        
        for kk=1:s2
            c(:,kk)=vec(( (1/sqrt(J)) *fft( repmat( (exp(1j*pi*( (0: (J-1))' ) ) ),1,D).*((sqrt(D)*ifft( (vec2mat( b(:,kk),D )).') ).') )).');
        end
        y=c*sqrt(J*D/length(pilot_positions));
        
        
    end
end