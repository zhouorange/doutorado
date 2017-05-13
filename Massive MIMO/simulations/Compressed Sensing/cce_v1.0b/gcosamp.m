function Sest=gcosamp(u,Phi,Phisel,groups,K,maxit)

%%% G-Cosamp algorithm 
%%% u: obs vector
%%% Phi : measurement matrix (function handle, such that Phi(u,1)=Phi*u and
%%%       Phi(u,2)=Phi.'*u)
%%% Phisel : function handle for selecting columns of Phi
%%% groups: group indices
%%% K : sparsity of Sest
%%% maxit : maximal number of iterations 
%%%
%%% Written by David Mary
%%%
%%% This script/program is released under the Commons Creative Licence
%%% with Attribution Non-commercial Share Alike (by-nc-sa)
%%% http://creativecommons.org/licenses/by-nc-sa/3.0/
%%% Disclaimer: the short answer, this script is for educational purpose only.
%%% More Disclaimer:  http://igorcarron.googlepages.com/disclaimer
%%%
%%% Modified by Daniel Eiwen, Nov 2010

% Initialization
N=length(Phi(u,2));
J=max(groups);
for jj=1:J
    gr{jj}=find(groups==jj);
end
Sest=zeros(N,1);
v=u;
t=1; T2=[];
phitr_gr=zeros(J,1);
while t < maxit 
    phitr=Phi(v,2);
    for jj=1:J
        phitr_gr(jj)=norm(phitr(gr{jj}));
    end
    [k,z]=sort(abs(phitr_gr),'descend');
    Omega_gr=z(1:2*K);
    T_gr=sort(union(Omega_gr,T2));
    T=[];
    for jj=1:length(T_gr)
        T=[T;gr{T_gr(jj)}];
    end
    Sest(T)=lsqr(@(x,mode)Phisel(x,T,mode),u,1e-15,1000);
    Sest_gr=zeros(J,1);
    for jj=1:J
        Sest_gr(jj)=norm(Sest(gr{jj}));
    end
    [k2,z2]=sort(Sest_gr,'descend');
    T2=z2(1:K);
    for jj=K+1:J
        Sest(gr{z2(jj)})=zeros(length(gr{z2(jj)}),1);
    end
    v=u-Phi(Sest,1);
    t=t+1;
end
