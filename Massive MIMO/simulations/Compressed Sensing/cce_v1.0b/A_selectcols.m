function y=A_selectcols(x,T,A,N,mat,mode)
%%%     Copyright: (c) Daniel Eiwen, Georg Tauböck, 2010-2012
if strcmp(mode,'notransp')
    xx=zeros(N,1);
    xx(T)=x;
    y=A(xx,mat,1);
else if strcmp(mode,'transp')
    y=A(x,mat,2);
    y=y(T);
    end
end