function [rc,ind,g]=rcosine(l,alpha)

%%%     Copyright: (c) Georg Tauböck, 2006-2010
%%%
%%%     Generate root raised cosine filter
%%%

g=-(1+alpha)/2+(1+alpha)/(l-1)*(0:(l-1));

rc=ones(l,1);
ind=intersect( find( (abs(g)>=(1-alpha)/2) ), find( (abs(g)<=(1+alpha)/2) ) );


rc(ind)=1/2*( 1 - sin( pi/alpha*(abs(g(ind))-1/2) ) );
