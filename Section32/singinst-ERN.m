%%% Single instance cavity for Erdos-Renyi networks in order to obtain
%%% the distribution of first eigenvalues for different values of N

% Initial conditions
kmin = 0;
kmax = 25;
c = 4;
J = 1/sqrt(c);

% Control of sizes
ensamble = 100;
ptslambda = 20;
experiments = 100;
b = 20;
ptn = 0;

Delta = 0.;
Ns = [256, 512, 1024, 2048, 4096, 8192, 16384];

% Allocation of arrays
lambdaexp = zeros(ensamble, length(Ns));
Lambda_test = linspace(5, 6, ptslambda);

for N = Ns
    ptn = ptn + 1;

    for ens = 1:ensamble

        meanrho = zeros(ptslambda, 1);

        % Degrees and number of fields
        degrees = randpoissarray(N, kmin, kmax, c);
        M =  sum(degrees);

        % Network
        Adj = kabashimagen(N, degrees);
        Adj_modified = BinomialMatrix(Adj, J, Delta);

        % Neighbours
        neighbours = (M+1)*ones(kmax-1, M+1);
        Js = zeros(kmax-1, M+1);
        Aux = ones(M+1, 1);
        Aux(M+1) = 0;

        ptlink = 0;
        for node = 1:N
            [n, ~] = find(Adj_modified(:, node));
            for l = 1:length(n)
                ptlink = ptlink + 1;
                neighbours(1:length(n)-1, ptlink)= n([1:l-1 l+1:end]);
                Js(1:length(n)-1, ptlink) = full(Adj_modified(node, n([1:l-1 l+1:end])))';
            end
        end

        % Allocation
        normH = zeros(experiments, ptslambda);
        T = zeros(experiments, ptslambda);
        rhoexp = zeros(b+1, ptslambda);

        indlambda = 0;

        % Non-backtracking extrapolation
        for lambda = Lambda_test
            indlambda = indlambda + 1;

            % Cavity fields
            A = rand(M+1, 1);
            H = randn(M+1, 1);

            % Iteration of fields
            for exp = 1:experiments
                A = lambda - sum(J^2./A(neighbours))';
                H = sum(Js.*H(neighbours)./A(neighbours))';
                normH(exp, indlambda) = norm(H(1:end-1));
                T(exp, indlambda) = mean(A(1:end-1));
            end

            bind = 0;
            for exp = experiments-b:experiments
                bind = bind + 1;
                if normH(exp, indlambda) ~= 0 && normH(exp-1, indlambda) ~= 0
                    rhoexp(bind, indlambda) = log(normH(exp, indlambda)/normH(exp-1, indlambda));
                end
            end

            meanrho(indlambda) = mean(rhoexp(:, indlambda));
        end

        p = polyfit(Lambda_test, meanrho', 1);
        lambdaexp(ens, ptn) = roots(p);

        if mod(ens, 10) == 0
            disp([ptn, ens]);
        end
    end
end

clearvars Adj aux count degrees Delta entries experiments exp ind indminus J K kaba l lambda M N n neighbours node pt ptn R2 s text Tri
