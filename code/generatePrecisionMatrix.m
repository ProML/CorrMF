function output = generatePrecisionMatrix(graph,density,nNode)
% data.sim <- bdgraph.sim(n = 100, p = 8, size = 10, vis = T)
%     (n = 1, p = 10, graph = "random", size = NULL, vis = FALSE)
if strcmp(graph,'random')
    nEdge = floor(density*nNode*(nNode-1)/2);
    G = zeros(nNode, nNode);
    sampleIndex = randperm(nNode*nNode, nEdge);
    G(sampleIndex) = 1;
    G = G + G.';
    G(G>0) = 1;
end
if strcmp(graph,'circle')
    G = zeros(nNode,nNode);
    K = eye(nNode);
    for i=1:(nNode - 1)
        K(i, i + 1) = 0.5;
        K(i + 1, i) = 0.5;
        K(1, nNode) = 0.4;
        K(nNode, 1) = 0.4;
        G = ceil(K);
    end
end
G(logical(eye(nNode)))= 0;
u = 0.1;
v = 0.3;
K = G * v;
[~,eigenValues] = eig(K);
K(logical(eye(nNode))) = abs(min(min(eigenValues))) + 0.1 + u;
sigma = corrcov(inv(K)); %cov2cor(solve(K));
K = inv(sigma);
% K = -K;
% K(logical(eye(nNode))) = abs(min(min(eigenValues))) + 0.1 + u;
isPositiveDefinite_K = all(eig(K) > 0)
isPositiveDefinite_sigma = all(eig(sigma) > 0)
output.precision = K;
output.covariance = sigma;
output.graph = G;

% Phi = wishrnd(Psi,nu);
% Phi = [
%   1.0,  -0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  -0.4;
%   -0.5,  1.0,  -0.5,  0.0,  0.0,  0.0,  0.0,  0.0;
%   0.0,  -0.5,  1.0,  -0.5,  0.0,  0.0,  0.0,  0.0;
%   0.0,  0.0,  -0.5,  1.0,  -0.5,  0.0,  0.0,  0.0;
%   0.0,  0.0,  0.0,  -0.5,  1.0,  -0.5,  0.0,  0.0;
%   0.0,  0.0,  0.0,  0.0,  -0.5,  1.0,  -0.5,  0.0;
%   0.0,  0.0,  0.0,  0.0,  0.0,  -0.5,  1.0,  -0.5;
%   -0.4,  0.0,  0.0,  0.0,  0.0,  0.0,  -0.5,  1.0];
% Phi = [
% 2.01E+00,	5.84E-01,	5.91E-01,	5.94E-01,	-1.13E-16,	6.03E-01	1.94E-16	-3.57E-17
% 5.84E-01,	1.50E+00,	7.25E-17,	7.22E-17,	5.21E-17,	8.15E-17	4.95E-01	9.55E-18
% 5.91E-01,	1.71E-16,	1.54E+00,	3.23E-17,	-7.66E-17,	5.28E-01	5.01E-01	-6.03E-17
% 5.94E-01,	2.15E-16,	7.68E-17,	1.55E+00,	5.50E-01,	1.47E-16	4.01E-17	7.36E-17
% 0.00E+00,	-5.59E-17,	-3.70E-17,	5.50E-01,	1.72E+00,	5.59E-01	-3.97E-17	4.83E-01
% 6.03E-01,	1.37E-16,	5.28E-01,	-3.88E-33,	5.59E-01,	1.60E+00	-2.23E-17	-2.65E-17
% 1.24E-16,	4.95E-01,	5.01E-01,	1.08E-16,	3.07E-17,	7.93E-17	1.44E+00	-2.53E-18
% -1.66E-17,	-1.56E-17,	-1.71E-17,	9.95E-17,	4.83E-01	-7.60E-35	-1.38E-18	1.19E+00
% ];
% Phi = [
% 1.66E+00,	5.73E-01,	1.50E-17,	5.86E-01,	5.78E-01,	-4.36E-17,	4.88E-17,	5.32E-01
% 5.73E-01,	1.66E+00,	-1.40E-17,	5.86E-01,	5.78E-01,	7.45E-17,	-3.49E-18,	5.32E-01
% -1.40E-17,	-2.80E-17,	1.17E+00,	2.80E-17,	-2.59E-17,	-4.36E-18,	4.42E-01,	-5.96E-18
% 5.86E-01,	5.86E-01,	4.23E-17,	1.73E+00,	-1.52E-16,	-2.29E-17,	5.38E-01,	-2.13E-17
% 5.78E-01,	5.78E-01,	-1.01E-17,	1.85E-17,	1.69E+00,	4.92E-01,	2.69E-17,	1.19E-16
% 1.23E-18,	-3.21E-17,	-2.96E-18,	1.81E-17,	4.92E-01,	1.20E+00,	1.18E-17,	3.30E-18
% -7.88E-17,	-1.57E-16,	4.42E-01,	5.38E-01,	-6.08E-17,	-1.77E-17,	1.40E+00,	-5.12E-17
% 5.32E-01,	5.32E-01,	1.17E-17,	2.13E-17,	-7.17E-17,	-3.14E-17,	5.34E-17,	1.43E+00];



