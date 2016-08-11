function [SDR_s,SIR_s,SAR_s,rec_test_s] = sep_with_framework(Tr_s,Tr_n,test_s,test_n,X,param,seeds)

k_s = [];
iter = 100;
while(isempty(k_s))
k_s = bin_search_dim(Tr_s,Tr_n,iter,param,seeds);     % dimension for source model
param.r_th = param.r_th - 0.2;
param.r_th = param.r_th*10;
param.r_th = round(param.r_th);
param.r_th = param.r_th/10;
end

d(1) = k_s;
D_s = nmf_bases(Tr_s,k_s,iter,seeds.b1,seeds.g1);    % final NMF dictionary for the source

k = 0;
disc_dim_vec = 20:5:120;

for disc_dim = disc_dim_vec(1):disc_dim_vec(2)-disc_dim_vec(1):disc_dim_vec(end)
    k = k+1;
    rng(seeds.b2);
    D_n = rand(size(Tr_n,1),disc_dim);
    rng(seeds.g2);
    C_nn = rand(disc_dim, size(Tr_n,2));
    for iter = 1:100
        D_n = update_B(D_n, C_nn, Tr_n, D_s, param.lambda);
        C_nn = update_G(D_n, C_nn, Tr_n);
    end
    
 D = [D_s, D_n];
 
 C_s = rand(size(D,2), size(Tr_s,2));
 C_n = rand(size(D,2), size(Tr_n,2));
 for i = 1:100
     C_s = update_G(D,C_s,Tr_s);
     C_n = update_G(D,C_n,Tr_n);
 end
 
ener1(1,k) = (norm(D_s*C_s(1:size(D_s,2),:),'fro')^2)/(size(Tr_s,2));
ener1(2,k) = (norm(D_n*C_s(size(D_s,2) + 1:end ,:),'fro')^2)/(size(Tr_s,2));
ener2(1,k) = (norm(D_n*C_n(size(D_s,2) + 1:end ,:),'fro')^2)/(size(Tr_n,2));
ener2(2,k) = (norm(D_s*C_n(1:size(D_s,2),:),'fro')^2)/(size(Tr_n,2));

end
r_ener1=ener1(1,:)./ener1(2,:);
r_ener2=ener2(1,:)./ener2(2,:);

ind = 1;
dim = [];
param.r1_min = 6; param.r2_max = 30;

while (ind <= numel(disc_dim_vec) && isempty(dim))
    if r_ener1(ind) > param.r1_min && r_ener2(ind) < param.r2_max
        ind = ind + 1;
    elseif ind == 1
        dim = disc_dim_vec(ind);
    else
        dim = disc_dim_vec(ind-1);
    end
end
if (isempty(dim))
    dim = disc_dim_vec(end);
end
d(2) = dim;

D_n = []; C_nn = [];
rng(seeds.b2);
D_n = rand(size(Tr_n,1),d(2));
rng(seeds.g2);
C_nn = rand(d(2), size(Tr_n,2));
for iter = 1:100
    D_n = update_B(D_n, C_nn, Tr_n, D_s, param.lambda);
    C_nn = update_G(D_n, C_nn, Tr_n);
end
 D = [D_s, D_n];
 
 rng(seeds.g);
 G = rand(size(D,2),size(X,2));
 for i = 1:100
     G = update_G(D,G,abs(X));
 end
 
 % initially estimated magnitude spectograms
S1_in=D_s*G(1:size(D_s,2),:);
S2_in=D_n*G(size(D_s,2)+1:end,:);

% spectral masks
H1=S1_in./(S1_in + S2_in);
H2=S2_in./(S1_in + S2_in);

S1_h=H1.*X;
S2_h=H2.*X;

[rec_test{1}, ~] = istft(S1_h,param.h,param.nfft,param.fs);
[rec_test{2}, ~] = istft(S2_h,param.h,param.nfft,param.fs);

test_s = (test_s(1:length(rec_test{1})))';
test_n = (test_n(1:length(rec_test{2})))';

se=[];
s=[];

se(1,:) = rec_test{1};
se(2,:) = rec_test{2};
s(1,:) = test_s;
s(2,:) = test_n;

[SDR,SIR,SAR,~] = bss_eval_sources(se,s);

% separation quality of source
SDR_s = SDR(1);
SIR_s = SIR(1);
SAR_s = SAR(1);
rec_test_s = rec_test{1};