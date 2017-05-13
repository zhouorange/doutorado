function y=A_optimized_ext(z,pil_pos,U,delay,dop,normalize_vec,J,mode)

%%%     Copyright: (c) Georg Tauböck, Daniel Eiwen, 2006-2012
%%%
%%%     Generates CS measurement matrix


lz=length(z)/J;

if ceil(lz)~=lz
    disp('Error in A_optimized_big!');
end    

if mode==1
    
    z=z.*normalize_vec(:);
    
    for jj=1:J
        a(:,jj)=vec((U*(((1/sqrt(delay))*fft((vec2mat((z(((jj-1)*lz+1):(jj*lz))),delay)).')).')).');
        y(:,jj)=a(pil_pos(:,jj),jj);
    end
    
    y=y(:);
    
else if mode==2
        
        a=zeros(lz,J);
	    a(:)=z;
        b=zeros(delay*dop,J);
        
        for jj=1:J
            b(pil_pos(:,jj),jj)=a(:,jj);
            y(:,jj)=vec((U'*((sqrt(delay)*ifft((vec2mat(b(:,jj),delay)).')).')).');
        end
        
        y=y.*normalize_vec;
        y=y(:);
    end
end