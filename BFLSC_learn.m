function [G] = BFLSC_learn(X, w, n_iter,lambda1,beta,gamma)
% Find a projection matrix 
% Input: 
%       X: N*1 LSM vector
%       w: binary code length
%       n_iter: number of iterations
%       lambda1, beta: control parameters
% Output:
%       G: d*K Projection matrix 


%--- Initialize
[nSmp, nFea] = size(X);
M = mean(X,1);
data = (X - repmat(M,nSmp,1));
[eigvec, eigval] = eig(data'*data);
[~,I] = sort(diag(eigval),'descend');
Go = eigvec(:,I(1:w));
A = [eye(w-1),zeros(w-1,1)] - [zeros(w-1,1),eye(w-1)];

%--- 2-stage method
opts.record = 0;
opts.mxitr  = 1000;
opts.xtol = 1e-5;
opts.gtol = 1e-5;
opts.ftol = 1e-8;

G = Go;
M = repmat(M,nSmp,1);
Xtmp = X'*ones(nSmp,1);
Q = (gamma*(-X'*X - 2*X'*M + M'*M) + beta*(Xtmp*Xtmp'))/nSmp;
tic;
iter=1;
for iter=1:n_iter
    CC(iter)=norm( G,2);
         kk(iter)=iter;
         iter=iter+1;
    % fix G. update B
    B = double(X*G >0) - 0.5;
    
    % fix B. update G.
    F1 = Q + lambda1*(X'*X)/nSmp;
    F2 = (-lambda1*2*X'*B - beta*X'*ones(nSmp,1)*ones(1,w))/nSmp; 
    F3 = A'*A;
    F4 = X'*X/nSmp;
    [G, ~]= OptStiefelGBB(G, @myfun,opts,F1,F2',F3,F4);
    
    
end
toc;


end







function [F, G] = myfun(G,M,N,P,O)
    F = trace(G'*M*G) + trace(N*G) + trace(G'*O*G*P*G'*O*G*P) - 2*trace(G'*O*G*P);
    G = 2*M*G + N' + 2*(O*G*P + O'*G*P')*(G'*O*G*P)' - 2*(O*G*P +O'*G*P');
end
