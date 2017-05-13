function y = A_Fourier_ext(z,pil_pos,delay,dop,J,mode)

%%%     Generates CS Measurement Matrix
%%%
%%%     Copyright: (c) Georg Tauböck, Daniel Eiwen, 2006-2012

lp=length(pil_pos(:,1));

lz=length(z)/J;
if ceil(lz)~=lz
    disp('Error in A_Fourier_ext!');
end    

if mode==1
    
    for jj=1:J
        a(:,jj)=vec((repmat((exp((-j)*pi*((0:(dop-1))'))),1,delay).*(ifft(((sqrt(dop)/sqrt(delay))*fft((vec2mat(z(((jj-1)*lz+1):(jj*lz)),delay)).')).'))).');
        y(:,jj)=a(pil_pos(:,jj),jj)*sqrt(dop*delay/lp);
    end
    y=y(:);
    
    
else if mode==2
        
        a=zeros(lz,J);
	 a(:)=z;
        b=zeros(delay*dop,J);
        
        for jj=1:J
            b(pil_pos(:,jj),jj)=a(1:lp,jj);
            c(:,jj)=vec(((1/sqrt(dop))*fft(repmat((exp(j*pi*((0:(dop-1))'))),1,delay).*((sqrt(delay)*ifft((vec2mat(b(:,jj),delay)).')).'))).');
        end
        y=c*sqrt(dop*delay/lp);
        y=y(:);
        
    end
end