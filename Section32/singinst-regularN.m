%%% Single instance cavity in order to obtain the distribution of first eigenvalues
%%% for regular graphs for different values of N.

% Initial conditions
K = 4;
J = 1;
Delta = 0;
Lambdath = 2*J*sqrt(K-1)+0.1;

ensamble = 100;
pts = 20;
experiments = 1000;
b = 200;
ptn = 0;

% Allocation of arrays
Ns = [256, 512, 1024, 2048, 4096, 8192, 16384];
Lambda_test = linspace(Lambdath, Lambdath+0.5, pts);

lambdaens = zeros(ensamble, length(Ns));

for N = Ns
    ptn = ptn + 1;

    % Degrees and number of fields
    degrees = K*ones(N, 1);
    M = sum(degrees);

    for ens = 1:ensamble

        % Network
        kaba = kabashimagen(N, degrees);

        % Entries
        Tri = triu(kaba, 1);
        entries = find(Tri);
        aux = rand(length(entries), 1);

        indminus = find(aux<=(1-Delta)/2);
        ind = entries(indminus);

        Tri(ind) = -1;
        Adj = Tri+Tri';

        % Neighbours and J's
        neighbours = [];
        Js = [];

        for node = 1:N
            [n, ~] = find(Adj(:, node));
            for l = 1:length(n)
                neighbours = [neighbours n([1:l-1 l+1:end])];
                Js = [Js full(Adj(node, n([1:l-1 l+1:end])))'];
            end
        end

        % Allocation
        normH = zeros(experiments, pts);
        rhoexp = zeros(b+1, pts);
        meanrho = zeros(pts, 1);

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
        p = polyfit(Lambda_test, meanrho', 2);
        rts = roots(p);

        lambdaens(ens, ptn) = rts(2);

        if mod(ens, 10) == 0
            disp([ptn, ens]);
        end
    end
end

clearvars Adj aux count degrees Delta entries experiments exp ind indminus J K kaba l lambda Lambda_test M N n neighbours node pt ptn R2 s text Tri
