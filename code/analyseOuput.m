% Calculate correlation from precision matrix (GMRF)
for i=1:size(model.meanPhi,1)
    for j=1:size(model.meanPhi,2)
        correlation(i,j) = -model.meanPhi(i,j)/sqrt(model.meanPhi(i,i)*model.meanPhi(j,j));
    end
end

corr = round(correlation*10)/10;
adj = corr;
% adj(adj~=0) = 1;
adj(abs(corr)>0.1) = 1;
adj(abs(corr)<=0.1) = 0;

nNode = size(adj,1);
nEdge = nNode*(nNode-1)/2;

adjacencyMatrix.true = output.graph;
adjacencyMatrix.true(logical(eye(size(output.graph,1)))) = 1;
error = adj-adjacencyMatrix.true;
FP = sum(sum(error(error>0))); % False positive
FN = -sum(sum(error(error<0))); % False negative
TP = sum(sum(adj & adjacencyMatrix.true));
TN = sum(sum(~(adj | adjacencyMatrix.true)));
precision = TP/(TP+FP);
recall = TP/(TP+FN);
F1 = 2*(precision*recall)/(precision+recall);
FPR = FP/(FP+TN);
density = sum(sum(triu(adj,1)))/nEdge;
error_acc = sum(sum(abs(adj-adjacencyMatrix.true)))/2;
accuracy = 1-(error_acc/nEdge);

%% Analyse matrix S
binaryS.pred = bin(model.meanS,0.5);
binaryS.true = bin(data.S,0.5);
error_S = sum(sum(abs(binaryS.pred-binaryS.true)))/(D*P);
