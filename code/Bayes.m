% Bayes:
% data: observed data
%  X: G*D gene expression-disease/drug continuous matrix
%  K: G*P gene-pathway binary matrix (prior knowledge)
% model: model parameters
%  tauEpsilon : random noise
%  B: G*P gene-pathway continuous matrix
%  tauB :
%  S: P*D pathway-disease/drug continuous matrix
%  Phi: P*P covariance matrix of S (pathway-pathway correlation)
%  hyperparameters
%    alphaEpsilon: for Gamma distribution
%    betaEpsilon: for Gamma distribution
%    alphaB: for Gamma distribution
%    betaB: for Gamma distribution
%    Eta: 2*1 matrix
%    nu: for Wishart distribution
%    Psi: P*P matrix for Wishart distribution
%    Pi: G*P matrix
% control: sampling paramters
%   iterNo: number of MCMC iterations
%   sampleNo: number of saved samples
%   burnNo: number of burn-in iterations
%   showIter: report every 'showIter' iterations
%   init: if true, the algorithm uses user-specified starting values
% post: posterior results
%
% Naruemon Pratanwanich - 2013
%
function [model] = Bayes(data,model0,control)
% Setting internal parameters
[G,D] = size(data.X);
[~,P] = size(data.K);
% Initialisation
model0.hyperparameters.Pi = initializePi();
current.value.Z = data.K;
if control.init
    current.value.tauEpsilon = model0.tauEpsilon;
    current.value.B = model0.B;
    current.value.tauB = model0.tauB;
    current.value.Phi = model0.Phi;
    current.value.S = model0.S;
else
    current.value.tauEpsilon = 1;
    current.value.B = randn(G,P);
    current.value.tauB = 1;
    current.value.S = randn(P,D);
    current.value.Phi = cov(current.value.S);
end
%% Gibbs sampling
for t=1:control.iterNo
    if mod(t,control.showIter) == 0
        t
    end
    % sample parameters
    [current.value.tauEpsilon] = sample_tauEpsilon();
    [current.value.tauB] = sample_tauB();
%     [current.value.tauB] = 1;
%     for i=1:G
%         for j=1:P
%             current.value.B(i,j) = sample_b_gp(i,j);
%         end
%     end
        for i=1:G
            current.value.B(i,data.K(i,:)==1) = sample_b_g(i);
            current.value.B(i,data.K(i,:)==0) = 0;

        end
    current.value.B = normalizeRow(current.value.B);
    
    current.value.Phi = sample_Phi();
%     current.value.Phi = solve_Phi();
%     current.value.Phi = eye(P);
    for i=1:D
        current.value.S(:,i) = sample_s_d(i);
    end
    %% Convergence diagnosis
    % value-iteration graphs for all parameters
    samples.value.tauEpsilon(t) = current.value.tauEpsilon;
    samples.value.tauB(t) = current.value.tauB;
    samples.value.Phi(t) = current.value.Phi(1,2);
    samples.value.S(t) = current.value.S(1,1);
    samples.value.B(t) = current.value.B(1,1);
    
   
    subplot(3,2,[1 2]);plot(1:t,samples.value.tauEpsilon);title('value tau epsilon');
    subplot(3,2,3);plot(1:t,samples.value.tauB);title('value tau B');
    subplot(3,2,4);plot(1:t,samples.value.B);title('a value in B');
    subplot(3,2,5);plot(1:t,samples.value.Phi);title('value Phi');
    subplot(3,2,6);plot(1:t,samples.value.S);title('a value in S');
    drawnow
    
    
    % Save samples
    if t>control.burnNo && mod(t-control.burnNo,control.sampleNo)==0
        sampleNo = (t-control.burnNo)/control.sampleNo;
        model.tauEpsilon(sampleNo) = current.value.tauEpsilon;
        model.Phi(:,:,sampleNo) = current.value.Phi;
        model.S(:,:,sampleNo) = current.value.S;
        model.B(:,:,sampleNo) = current.value.B;
        model.tauB(sampleNo) = current.value.tauB;
        model.Z(:,:,sampleNo) = current.value.Z;
        %         disp('save');
    end
    
end
disp('Done Gibbs Sampling');

%% Output Manipulation
model.meanS = mean(model.S,3);
model.meanPhi = mean(model.Phi,3);
model.meanB = mean(model.B,3);
model.meanTauB = mean(model.tauB);
model.meanTauEpsilon = mean(model.tauEpsilon);
model.Pi = model0.hyperparameters.Pi;
%% Functions for sampling
    function [tauEpsilon] = sample_tauEpsilon()
        alpha = model0.hyperparameters.alphaEpsilon + (G*D)/2;
        beta = model0.hyperparameters.betaEpsilon + 0.5*sum(sum((data.X-(current.value.B*current.value.S)).^2));
        tauEpsilon = gamrnd(alpha,1/beta);
    end

    function [b_g] = sample_b_g(g)
        subMatrix_S = current.value.S((data.K(g,:)==1),:);
        sub_P = size(subMatrix_S,1);
%         subMatrix_S = current.value.S;
%         sub_P = P;
        Tau = current.value.tauB*ones(sub_P,sub_P) + current.value.tauEpsilon*(subMatrix_S*subMatrix_S');
        mu = Tau\eye(sub_P)*current.value.tauEpsilon*(subMatrix_S*data.X(g,:)');
        sigma = Tau\eye(sub_P);
%         sigma = round(sigma*10000)/10000; %bad
        sigma = 0.5.*(sigma+sigma'); % make sure for symmetric 
       
        if sub_P>1 
        b_g = mvnrnd(mu,sigma);
        elseif sub_P==1
            b_g = normrnd(mu,sigma);
        else
            b_g = 0;
        end
    end
    function [b_gp] = sample_b_gp(g,p)
        %         sumRest = 0;
        %         for k=1:P
        %             if k~=p
        %                 sumRest = sumRest + current.value.B(g,k)*current.value.S(k,:);
        %             end
        %         end
        %         x_gp = data.X(g,:) - sumRest;
        current.value.X = current.value.B*current.value.S;
        x_gp = data.X(g,:) - current.value.X(g,:) + current.value.B(g,p)*current.value.S(p,:);
        
        s_p = current.value.S(p,:)'; 
        x_gp = x_gp'; 
        
        BF_gp = calculateBF(x_gp,s_p);
        Pi_gp = model0.hyperparameters.Pi(g,p)/(BF_gp*(1-model0.hyperparameters.Pi(g,p))+model0.hyperparameters.Pi(g,p));
        r = rand(1);
        if r>Pi_gp % draw from the Dirac delta function
            %         if data.K(g,p)==0
            b_gp = 0;
            current.value.Z(g,p) = 0;
        else % draw from f(b)
            tau = current.value.tauEpsilon*(s_p'*s_p) + current.value.tauB;
            mu = (current.value.tauEpsilon/tau)*(s_p'*x_gp);
            sigma = 1/tau;
            b_gp = normrnd(mu,sigma);
            current.value.Z(g,p) = 1;
        end
    end
    function [tauB] = sample_tauB()
        alpha = model0.hyperparameters.alphaB + 0.5*sum(sum(current.value.Z));
        beta = model0.hyperparameters.betaB + 0.5*((norm(current.value.B,'fro'))^2);
        tauB = gamrnd(alpha,1/beta);
    end

    function [Phi] = sample_Phi()
        sumS = zeros(P,P);
        for d=1:D
            sumS = sumS + current.value.S(:,d)*current.value.S(:,d)';
        end
        nu = model0.hyperparameters.nu + D;
        Psi = inv(model0.hyperparameters.Psi) + sumS;
        Phi = wishrnd(inv(Psi),nu);
    end
    function Phi = solve_Phi()
        empiricalSigma = cov(current.value.S');
        [Phi ~] = graphicalLasso(empiricalSigma, model0.rho, model0.maxIterForPhi, model0.tolThreshold);
    end
    function [s_d] = sample_s_d(d)
        Phi = current.value.Phi + current.value.tauEpsilon*(current.value.B'*current.value.B);
        mu = (Phi\eye(P))*current.value.tauEpsilon*(current.value.B'*data.X(:,d));
        sigma = inv(Phi);
        s_d = mvnrnd(mu,sigma)';
    end


%% Auxilary functions
    function [Pi] = initializePi()
        model0.hyperparameters.Pi = zeros(G,P); % Calculate Pi Matrix depending on the matrix K
        for g=1:G
            for p=1:P
                if data.K(g,p)==1
                    model0.hyperparameters.Pi(g,p) = model0.hyperparameters.Eta(2);
                else
                    model0.hyperparameters.Pi(g,p) = 1-model0.hyperparameters.Eta(1);
                end
            end
        end
        Pi = model0.hyperparameters.Pi;
    end
    function [BF_gp] = calculateBF(x_gp,s_p)
        mu = zeros(D,1);
        sigma1 = (1/current.value.tauEpsilon)*eye(D);
        sigma2 = (1/current.value.tauB)*(s_p*s_p')+ sigma1;
        probModel1 = mvnpdf(x_gp,mu,sigma1);
        probModel2 = mvnpdf(x_gp,mu,sigma2);
        if probModel2 ==0 && probModel1 == 0
            BF_gp = 1;
            disp('OOPS !!');
        else
            BF_gp = probModel1/probModel2;
        end
    end

end