%%% modgdcssomp_mode.M -  Performs the Modified Group-DCS-SOMP Algorithm for partially 
%%%                       known support, which recovers a jointly groupsparse ensemble 
%%%                       of signals from far fewer measurements than usually necessary, 
%%%                       where the support is partially known (complex-valued). 
%%%                       The measurement matrices can be given as function-handle
%%%
%%% Input          :   Y   =  a matrix containing the measurements as columns
%%%                    G   =  either stack of measurement matrices such that
%%%                           G(:,:,mat) is the mat-th matrix, or function
%%%                           handle such that G(x,mat,1) = G_mat*x and 
%%%                           G(x,mat,2) = (G_mat)'*x, where G_mat is the
%%%                           mat-th matrix
%%%                           ATTENTION: columns must have norm 1
%%%                    groups= vector of group numbers
%%%                    maxig= maximal number of (group) iterations
%%%                    supp=  known support (vector of group numbers)
%%%                    e   =  target error power level for convergence 
%%%
%%% Output         :   Xrec=  the recovered signals, as columns of a matrix
%%%                           columns of a matrix
%%%                    grind= indices of chosen groups
%%%
%%% Usage          :   [Xrec,grind] = modgdcssomp(Y,G,group,maxig,supp,e) ;
%%%
%%% COPYRIGHT: (c) Daniel Eiwen, 2010-2012

function [Xrec,grind]=modgdcssomp(Y,G,group,maxig,supp,e)

if nargin<3
    help modgdcssomp;
    Xrec=[];
    return;
end

[M,J]=size(Y);
N=length(group);
K=max(group);
lg=N/K;

if isa(G,'function_handle')
    I=eye(N);
    Gcol=@(c,mat) G(I(:,c),mat,1);
end

if ~isa(G,'function_handle')
    GG=G;
    clear G;
    G=@(x,mat,mode) GG(:,:,mat)'*x;
    Gcol=@(c,mat) GG(:,c,mat);
end



if nargin<4
    maxig=floor(M/lg);
end
if nargin<5
    e=0.001;
end

maxi=lg*maxig;

for jj=1:K
    gr{jj}=find(group==jj);
end


%%% Initialization
z=1;
Gam=zeros(M,maxi,J);
R=Y;
grind=zeros(1,maxi);
B=zeros(maxi,J);
rr=zeros(maxi,maxi,J);
nres=zeros(maxi,J);
c=0;
L=[];

% %%% Check if all cloumns have unit l2-norm (takes quite some time)
% for kk=1:J
%    for ll=1:N
%        if abs(norm(Gcol(ll,kk))-1)>1000*eps
%             disp('Columns of the measurement matrix have to be unit length!');
%             return;
%        end
%    end
% end

%%% use known support first

lsupp=length(supp);
for nl=1:lsupp
    c=c+1;
    IO = gr{supp(nl)};
    L = [L' IO']';
    grind(c)=supp(nl);
end

%%% QR-decomposition
for jj=1:J
    for kk=1:lsupp*lg
        g=Gcol(L(kk),jj);
        if kk>1
            a=g;
            for tt=1:kk-1
                gtj=Gam(:,tt,jj);
                rr(tt,kk,jj)=gtj'*a;
                g=g-rr(tt,kk,jj)*gtj;
            end
        end
        rr(kk,kk,jj)=norm(g);
        g=g/rr(kk,kk,jj);
        Gam(:,kk,jj)=g;
        B(kk,jj)=g'*R(:,jj);
        R(:,jj)=R(:,jj)-B(kk,jj)*g;
    end    
end


%%% Main Loop

while z==1 & c<maxig
    c=c+1;

    %%% Finding the next atom
    H=zeros(N,J);
    for jj=1:J
        H(:,jj)=G(R(:,jj),jj,2);
    end
    HH=zeros(K,1);
    for mm=1:K
        HH(mm)=norm(H(gr{mm},:),'fro');
    end
    nl=find(HH==max(HH),1);
    I = gr{nl};
    L = [L' I']';
    grind(c)=nl;

    %%% QR-decomposition
    for jj=1:J
        kk=(c-1)*lg;
        for mm=1:lg
        kk=kk+1;
        g=Gcol(L(kk),jj);
        if kk>1
            a=g;
            for tt=1:kk-1
                gtj=Gam(:,tt,jj);
                rr(tt,kk,jj)=gtj'*a;
                g=g-rr(tt,kk,jj)*gtj;
            end
        end
        rr(kk,kk,jj)=norm(g);
        g=g/rr(kk,kk,jj);
        Gam(:,kk,jj)=g;
        B(kk,jj)=g'*R(:,jj);
        R(:,jj)=R(:,jj)-B(kk,jj)*g;
        end

        %%% halting criterion using the norm of the residual
        nres(c,jj)=norm(R(:,jj));
        if nres(c,jj)<=e*norm(Y(:,jj))
            z=0;
        end
        
    end
end
B(kk+1:maxi,:)=[];
grind(c+1:maxi)=[];

Xrec=zeros(N,J);
for jj=1:J
    Xrec(L,jj)=rr(1:kk,1:kk,jj)\B(1:kk,jj);
end
