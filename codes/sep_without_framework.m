function [SDR,SIR,SAR,rec_test] = sep_without_framework(Tr,test,X,param,seeds)

% initialise basis and gain matrices
rng(seeds.b1);
B{1} = rand(size(Tr{1},1),128);
rng(seeds.b2);
B{2} = rand(size(Tr{2},1),128);
rng(seeds.g1)
G{1} = rand(128,size(Tr{1},2));
rng(seeds.g2)
G{2} = rand(128,size(Tr{2},2));

% update B,G
for iter=1:100
    B{1} = update_B(B{1},G{1},Tr{1},B{2},param.lambda);
    B{2} = update_B(B{2},G{2},Tr{2},B{1},param.lambda);
    G{1} = update_G(B{1},G{1},Tr{1});
    G{2} = update_G(B{2},G{2},Tr{2});   
end

% concatenated dictionaries
B_final = [B{1},B{2}];

% finding gain matrix for mixed signal
rng(seeds.g);
G = rand(size(B_final,2),size(X,2));

for i=1:100
    G = update_G(B_final,G,abs(X));
end

% initially estimated magnitude spectograms
S1_in = B{1}*G(1:size(B{1},2),:);
S2_in = B{2}*G(size(B{1},2)+1:end,:);

% spectral masks
H1=S1_in./(S1_in + S2_in);
H2=S2_in./(S1_in + S2_in);

S1_h=H1.*X;
S2_h=H2.*X;

[rec_test{1}, ~] = istft(S1_h,param.h,param.nfft,param.fs);
[rec_test{2}, ~] = istft(S2_h,param.h,param.nfft,param.fs);

test{1}=(test{1}(1:length(rec_test{1})))';
test{2}=(test{2}(1:length(rec_test{2})))';

se=[];
s=[];

se(1,:) = rec_test{1};
se(2,:) = rec_test{2};
s(1,:) = test{1};
s(2,:) = test{2};

[SDR,SIR,SAR,~] = bss_eval_sources(se,s);
