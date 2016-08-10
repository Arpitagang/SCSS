function [B] = nmf_bases(T,dim,iter,s_b1,s_g1)

rng(s_b1);
B = rand(size(T,1),dim);
rng(s_g1);
G = rand(dim,size(T,2));
Bp = rand(size(B));
for j = 1:iter
    B = update_B(B,G,T,Bp,0);
    G = update_G(B,G,T);
end

