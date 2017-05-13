function W = optimize_basis_jointgroup(gamma,g,J,L,K,N,reldop_max,reldop_min,Channels,scatterer_shift,blocksize)

%%% Copyright: (c) Georg Tauböck, 2006-2010
%%% Construct Optimized Basis Expansion Model which enforces joint sparsity
%%%
%%% Channels: number of channels (eg. for a 2x2 system = 4, 3x3 = 9, ...)
%%% scatterer_shift: number of taps the channels are shifted at the most
%%%                    (i.e. the scattering function of the 2nd, 3rd, 4th., ... channel is shifted at most +/- scatterer_shifts 
%%%                    in regard to the doppler of the 1st channel)
%%% blocksize: blocksize in discrete i-direction (exponent of 2)


Td=N*L;

k_max=ceil(reldop_max*Td/K);
k_min=floor(reldop_min*Td/K);

ScattererNumber=500;    %%can be varied

Amb_pre=repmat(conj(ambiguity_function(gamma,g,1,Td)),1,ScattererNumber);

I=eye(J);

F=(1/sqrt(J))*fft(I);

W=F;

V=I;

old_sparsity=Inf;

scatterers=zeros(Channels,ScattererNumber);
scatterers(1,:)=rand(1,ScattererNumber)*2*k_max-k_max*ones(1,ScattererNumber);
if Channels>1
    scatterers(2:Channels,:)=repmat(scatterers(1,:),Channels-1,1)+rand(Channels-1,ScattererNumber)*2*scatterer_shift-scatterer_shift*ones(Channels-1,ScattererNumber);
end    
H_matr3=[];

for chan=1:Channels
    indices=(0:(Td-1))'*scatterers(chan,:);
    H_matr1=exp(1j*2*pi*indices/Td);
    Spread_matr=fft(H_matr1);
    F_matr=Spread_matr.*Amb_pre;
    H_matr2=ifft(F_matr);
    H_matr3=[H_matr3 H_matr2(1:N:Td,:)];
end
H_matr3=H_matr3(1:J,:);

ell1_sparsity=norm(vec(W*H_matr3),1)
ell1_2_sparsity=sum(sqrt(sum((abs(reshape((reshape(W*H_matr3,J*ScattererNumber,Channels)).',blocksize*Channels,J*ScattererNumber/blocksize)).^2),1)))

ell1_sparsity_fourier=ell1_sparsity;
ell1_2_sparsity_fourier=ell1_2_sparsity;

schritt=1;

while schritt>2^(-20)
while ell1_2_sparsity < old_sparsity        
    
    W=V*W;    
    M=W*H_matr3;
    schritt_power=log2(schritt)

    cvx_solver sdpt3

    cvx_begin
        variable A(J,J) hermitian; 
        minimize(sum(norms(reshape((reshape((I+1j*2*pi*A)*M,J*ScattererNumber,Channels)).',blocksize*Channels,J*ScattererNumber/blocksize),2,1)));
        subject to
        max(max(abs(A))) <= schritt 
    cvx_end

    V=expm(1j*2*pi*A);
    
    clear scatterers indices H_matr1 Spread_matr F_matr H_matr2 H_matr3
    scatterers=zeros(Channels,ScattererNumber);
    scatterers(1,:)=rand(1,ScattererNumber)*2*k_max-k_max*ones(1,ScattererNumber);
    if Channels>1
        scatterers(2:Channels,:)=repmat(scatterers(1,:),Channels-1,1)+rand(Channels-1,ScattererNumber)*2*scatterer_shift-scatterer_shift*ones(Channels-1,ScattererNumber);
    end
    H_matr3=[];

    for chan=1:Channels
        indices=(0:(Td-1))'*scatterers(chan,:);
        H_matr1=exp(1j*2*pi*indices/Td);
        Spread_matr=fft(H_matr1);
        F_matr=Spread_matr.*Amb_pre;
        H_matr2=ifft(F_matr);
        H_matr3=[H_matr3 H_matr2(1:N:Td,:)];
    end
    H_matr3=H_matr3(1:J,:);

    old_sparsity2=old_sparsity;
    old_sparsity=ell1_2_sparsity
    ell1_sparsity=norm(vec(V*W*H_matr3),1)
    ell1_2_sparsity=sum(sqrt(sum((abs(reshape((reshape(V*W*H_matr3,J*ScattererNumber,Channels)).',blocksize*Channels,J*ScattererNumber/blocksize)).^2),1)))
end
schritt=schritt/2;
V=I;
ell1_2_sparsity=old_sparsity;
old_sparsity=old_sparsity2;
end


ell1_sparsity_fourier
ell1_sparsity_final=norm(vec(V*W*H_matr3),1)
ell1_2_sparsity_fourier
ell1_2_sparsity_final=sum(sqrt(sum((abs(reshape((reshape(V*W*H_matr3,J*ScattererNumber,Channels)).',blocksize*Channels,J*ScattererNumber/blocksize)).^2),1)))
