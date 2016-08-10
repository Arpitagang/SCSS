function B=update_B(B,G,S,Bp,lambda)

B=B.*(((S./(B*G))*G')./(ones(size(S))*G'+lambda*Bp*ones(size(Bp,2),size(B,2))));

for i = 1:size(B,2)
        B(:,i)=B(:,i)./norm(B(:,i));        
end