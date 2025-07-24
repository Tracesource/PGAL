function [A,Z,UU,iter,obj,obj_all] = PGAL(X,k,numanchor,lambda)
% m      : the number of anchor. the size of Z is m*n.
% X      : n*di

%% initialize
maxIter = 50 ; % the number of iterations
IterMax = 50;

m = numanchor;
numclass = k;
numview = length(X);
numsample = size(X{1},1);

%Calculate Mp
XX = [];
for p = 1 : numview
    X{p} = mapstd(X{p}',0,1);
    X{p} = X{p}';
    XX = [XX;X{p}'];
end
rand('twister',12);
[~,An] = litekmeans(XX',m, 'MaxIter',100,'Replicates',10);
diter = 1;
for p = 1:numview
    M{p} = An(:,diter:diter+size(X{p},2)-1);
    diter = diter+size(X{p},2);
%Initialize Wp
    W{p} = eye(m);
end

%Initialize gamma
gamma = ones(1,numview)./numview;

%Initialize Z
Z = computeIniGraph(XX,An',numclass);
Z = Z';
Z_sqr = Z.*Z;

flag = 1;
iter = 0;
obj = [];
obj_all = [];
%%
while flag
    iter = iter + 1;

    %% optimize Ap
    parfor p=1:numview
        MW = M{p}'*W{p};
        for jj=1:m
            beta = 0;
            c = zeros(1,size(X{p},2));
            for ii=1:numsample
                c = c+gamma(p)^2*Z_sqr(ii,jj)*X{p}(ii,:);
                beta = beta+gamma(p)^2*Z_sqr(ii,jj);
            end
            c = c+lambda*MW(:,jj)';
            beta = beta+lambda;
            A{p}(jj,:) = c'./beta;
        end
    end
    
    %% optimize Wp
    parfor p = 1:numview
        pp = A{p}*M{p}';
        [Unew,~,Vnew] = svd(pp','econ');
        W{p} = Unew*Vnew';
    end
    
    %% optimize gamma
    parfor p = 1:numview
        sumXA = 0;
        XA{p} = zeros(numsample,m);
        for ii = 1:numsample
            for jj = 1:m
                XA{p}(ii,jj) = norm(X{p}(ii,:)-A{p}(jj,:),2)^2;
                sumXA = sumXA + XA{p}(ii,jj)*Z_sqr(ii,jj);
            end
        end
        gamma(p) = 1/(sumXA);
    end
    gamma = gamma./sum(gamma);
    
    %% optimize Z
    parfor ii=1:numsample
        b = 0;
        for p = 1:numview
            b = b+gamma(p)^2*XA{p}(ii,:);
        end
        [Z(ii,:),~] = EProjSimplex_new_ZJP_V2(b,0);
    end
    Z_sqr = Z.*Z;

    compute_obj
    
    if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-4 || iter>maxIter || obj(iter) < 1e-10)
        [UU,~,~] = mySVD(Z,numclass);
        flag = 0;
    end
end
         
         
    
