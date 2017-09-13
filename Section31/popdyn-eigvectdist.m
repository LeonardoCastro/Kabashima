%%% Obtaining the distribution of entries of the first eigenvector for random regular graphs.
%%% First, q(A, H) is obtained.
%%% Then, Q(A, H) is computed.

% Initial conditions
K = 4;
Np = 10e6;
Nn = Np/K;
J = 1;
experiments = 50;

Delta = 0.9;
Delta_c = 1/sqrt(K-1);

% Eigenvalue
if Delta > Delta_c
  lambda = (K-1)*Delta+1/Delta;
else
  lambda = 2*J*sqrt(K-1);
end

% Thermodynamic value of cavity field A
Ath = (lambda + sqrt(lambda^2-4*(K-1)))/2;
itermax = 0;

% Cavity fields
H1 = rand(Np, 1); % Initial population
H2 = zeros(Nn, 1); % Eigenvector entries population

% neighbours and Js
Indices = randperm(Np);
neighbours = randomneighbours(Indices, K-1, Np);

Js = rand(K-1, Np);
ind = find(Js>(1-Delta)/2);
Js(ind) = J;
Js(Js~=1) = -J;

% Calculation of q(H)
for exp = 1:experiments
    H1 = sum(Js.*H1(neighbours)/Ath)';
end

Indices2 = randperm(Nn);
neighbours2 = randomneighbours(Indices2, K-1, Np);

% Calculation of Q(H)
for ind = 1:Nn
  H2 = sum(Js(:,1:Nn).*H1(neighbours2)/Ath)';
end

H2 = H2/std(H2);

clearvars Ath Delta Delta_c exp experiments H1 J Js ind Indices Indices2 itermax K lambda Nn Np neighbours neighbours2
