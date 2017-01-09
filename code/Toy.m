G = 150; 
D = 200;
P = 5; 
%% Setting model hyperparameter
model0.hyperparameters.alphaEpsilon = 1;
model0.hyperparameters.betaEpsilon = 1;
model0.hyperparameters.alphaB = 1;
model0.hyperparameters.betaB = 1;
model0.hyperparameters.Eta = [1,1]';
model0.hyperparameters.nu = P;
model0.hyperparameters.Psi = eye(P);

%% Generating toy data
alpha = model0.hyperparameters.alphaEpsilon;
beta = model0.hyperparameters.betaEpsilon;
tauEpsilon = gamrnd(alpha,beta);
tauEpsilon = 1;
sigma = 1/tauEpsilon;
epsilon = zeros(G,D);
epsilon = normrnd(0,sigma,G,D);

alpha = model0.hyperparameters.alphaB;
beta = model0.hyperparameters.betaB;
% tauB = gamrnd(alpha,beta);
tauB = 1;
mu = 0;
sigma = 1/tauB;
B = normrnd(mu,sigma,G,P);
K = ones(G,P);
randIdx = randperm(G*P);
missIdx = randIdx(1:round(G*P*0.95));
K(missIdx) = 0;
B(missIdx) = 0;
covBeforeNorm = cov(B(B~=0))
B = normalizeRow(B);
covAfterNorm = cov(B(B~=0))
S = zeros(P,D);
Psi = model0.hyperparameters.Psi;
nu = model0.hyperparameters.nu; 

nNode = P;
graph = 'random';
density = 0.5;
output = generatePrecisionMatrix(graph,density,nNode);
Phi = output.precision;

mu = zeros(P,1);
sigma = inv(Phi);
S = mvnrnd(mu,sigma,D);
S = S';
data.B = B;
data.tauB = tauB;
data.S = S;
data.Phi = Phi;
data.K = K;
data.tauEpsilon = tauEpsilon;
data.epsilon = epsilon;
data.X = data.B*data.S + data.epsilon;


%% Initializing model parameters
model0.tauEpsilon = gamrnd(model0.hyperparameters.alphaEpsilon,model0.hyperparameters.betaEpsilon);
model0.B = randn(G,P);
model0.tauB = gamrnd(model0.hyperparameters.alphaB,model0.hyperparameters.betaB);
model0.S = randn(P,D);
model0.Phi = wishrnd(model0.hyperparameters.Psi,model0.hyperparameters.nu);

model0.rho = 0.01;
model0.maxIterForPhi = 1000;
model0.tolThreshold = 0.0001;
%% Setting model controllers
control.iterNo = 5000;
control.sampleNo = 10;
control.burnNo = 2500;
control.showIter = 1;
control.init = true;
%% Running Gibbs sampling
tic;
[model] = Bayes(data,model0,control);
timeGibbs = toc;
