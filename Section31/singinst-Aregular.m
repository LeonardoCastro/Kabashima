%%% Single instance cavity in order to obtrain the distribution of the cavity field array A
%%% for random regular graphs and different values of N

% Initial conditions
J = 1;
K = 4;

Delta = 0.5;
Deltac = 1/sqrt(K-1);

experiments= 1000;

% Allocation of arrays
Ns = [256, 512, 1024, 1024*2, 1024*4, 1024*8, 1024*16];
As = zeros(K*Ns(end), length(Ns));

% Eigenvalue
if Delta < Deltac
    lambda = 2*J*sqrt(K-1);
else
    lambda = J*(Delta*(K-1)+1/Delta);
end

ptn = 0;
for N = Ns
    ptn = ptn + 1;
    degs = K*ones(N, 1);
    M = sum(degs);

    % Network
    kaba = kabashimagen(N, degs);

    % Entries
    Tri = triu(kaba, 1);
    entries = find(Tri);
    aux = rand(length(entries), 1);

    indminus = find(aux<=(1-Delta)/2);
    ind = entries(indminus);

    Tri(ind) = -1;
    Adj = Tri+Tri';

    % neighbours and entries
    neighbours = [];
    Js = [];

    for node = 1:N
        [n, ~] = find(Adj(:, node));
        for l = 1:length(n)
            neighbours = [neighbours n([1:l-1 l+1:end])];
            Js = [Js full(Adj(node, n([1:l-1 l+1:end])))'];
        end
    end

    % Cavity fields
    A = rand(M, 1);

    for exp = 1:experiments
        A = lambda - sum(1./A(neighbours))';
    end

    As(:, ptn) = A;
end

%%% Watch out with the zero entries in As

clearvars A Adj aux Delta Deltac Deltas degs entires ens ensamble exp experiments ind indminus J Js K kaba l lambda M meanA_local neighbours N n node pt ptd ptn pts sigmasquareA_local Tri
