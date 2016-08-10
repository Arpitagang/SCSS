function k_s = bin_search_dim(Tr,iter,param,seeds)

d = param.min_dim + round((param.max_dim-param.min_dim)/2);
k_s = [];

while((isempty(k_s)) && (param.min_dim <= param.max_dim))
    
B1 = nmf_bases(Tr{1},d,iter,seeds.b1,seeds.g1);

G11 = rand(size(B1,2),size(Tr{1},2));
G21 = rand(size(B1,2),size(Tr{2},2));
for i=1:iter
    G11 = update_G(B1,G11,Tr{1});
    G21 = update_G(B1,G21,Tr{2});
end

err(1)=mean(sqrt(sum((Tr{2}-B1*G21).^2,1)));
err(2)=mean(sqrt(sum((Tr{1}-B1*G11).^2,1)));

r = err(1)/err(2);
r = r*10;
r = round(r);
r = r/10;

if(r == param.r_th)
    k_s = d;
elseif(r > param.r_th)
    param.max_dim = d-1;
    d = param.min_dim + round((param.max_dim - param.min_dim)/2);
elseif(r < param.r_th)
    param.min_dim = d+1;
    d = param.min_dim + round((param.max_dim - param.min_dim)/2);
end
end
    
