%Marcos Bolanos
%November 2011
%Compressive Sensing Example

%This very simple example of L1 minimization is reproduced for
%implementation on matlab. The original example was posted on Rip's Applied
%Mathematics Blog on March 28, 2011 entitled "Compressed Sensing: the L1
%norm finds sparse solutions". 

%One needs to download the L1-MAGIC package in order to perform the l1
%minimization on matlab.

%This example was very good for illustrating how L1 minimization can
%identify a sparse vector. Here x is the sparse vector. A is the kxN
%incoherent matrix and B are the coefficients. The example shows how we can
%find the original x. xp should be approximately equal to x.

s=RandStream('mt19937ar');
RandStream.setDefaultStream(s);
reset(s);


 x=[0, 0, 0.319184, 0, 1.65857, 0, 0, 0, -1.0439, 0]'; %original sparse vector
 A=randn(s,10,5);
 B=A'*x;
 xp = l1eq_pd(x, A', 1, B)     %l1 minimization using L1-MAGIC
 