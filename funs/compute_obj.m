term1 = 0;
term2 = 0;
for p = 1:numview
    for ii = 1:numsample
        for jj = 1:m
            term1 = term1+ gamma(p)^2*norm(X{p}(ii,:)-A{p}(jj,:),2)^2*Z_sqr(ii,jj);
        end
    end
    term2 = term2 + norm(W{p}*A{p}-M{p},'fro')^2;
end

obj(iter) = term1+lambda*term2;
obj_all(end+1) = obj(iter);


