function [Y] = normalizeRow(X)

for r = 1:size(X,1)
    normCons = sum(abs(X),2);
    for c = 1:size(X,2)
        if X(r,c)~=0
            Y(r,c) = X(r,c)./normCons(r);
        end
    end
end