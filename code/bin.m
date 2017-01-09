function binary = bin(X,threshold)
ab = abs(X);
mean_ab = mean(ab);
for i=1:length(mean_ab)
binary(:,i) = 1./(1+exp(-ab(:,i)+mean_ab(i)));
end
binary(binary>threshold)=1;
binary(binary<=threshold)=0;