function G = update_G(B,G,S)

G = G.*((B'*(S./(B*G)))./(B'*ones(size(S))));