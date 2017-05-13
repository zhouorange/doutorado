function y=A_optimized(z,pil_pos,U,D,J,normalize_vec,mat,mode)

%%%     Copyright: (c) Georg Tauböck, 2006-2010
%%%
%%%     Generates CS measurement matrix


pilot_positions=pil_pos(:,mat);
s2=size(z,2);

if mode==1
    
    z=z.*repmat(normalize_vec(:,mat),1,s2);
    a=zeros(J*D,s2);
    
    for kk=1:s2
        a(:,kk)=vec((U*(((1/sqrt(D))*fft((vec2mat(z(:,kk),D)).')).')).');
    end
    
    y=a(pilot_positions,:);
    
else if mode==2
        a=z(1:(length(pilot_positions)),:);
        b=zeros(D*J,s2);
        b(pilot_positions,:)=a;
        y=zeros(J*D,s2);
        
        for kk=1:s2
            y(:,kk)=vec((U'*((sqrt(D)*ifft((vec2mat(b(:,kk),D)).')).')).');
        end
        
        y=y.*repmat(normalize_vec(:,mat),1,s2);
    end
end