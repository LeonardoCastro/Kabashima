% Code to obtain the distribution of first eigenvalues for random regular graph
% for different values of Delta and for a single system size.

% Initial conditions
N = 2048; % number of nodes
K = 4; % degree
J = 1; % value of the entries
degrees = K*ones(N, 1); % degree vector
M = sum(degrees); % twice number of links

ensamble = 100; % number of eigenvectors for each value of Delta
pts = 20; % number of values of Delta
experiments = 900; % number of iterations
b = 200; % number of rhos taken

Deltas = linspace(0. ,1., pts);
lambda_delta = zeros(pts, 1);
pt = 0;

for Delta = Deltas
    pt = pt + 1;

    % Depending of the value of Delta a different eigenvalue is taken
    if Delta < 1/sqrt(K-1)
        Lambdath = 2*J*sqrt(K-1);
    else
        Lambdath = J*((K-1)*Delta+1/Delta);
    end

    % Define interval for values of lambda
    lambdaexp = zeros(ensamble, 1);
    Lambda_test = linspace(Lambdath, Lambdath+1, pts);

    for ens = 1:ensamble

        % Allocation
        normH = zeros(experiments, 1);
        Ps = zeros(2, pts);
        count = 0;
        R2s = zeros(pts, 1);
        indlambda = 0;

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

        % For each value of lambda in the allowed interval
        for lambda = Lambda_test
            indlambda = indlambda + 1;

            % Cavity fields
            A = rand(M, 1);
            H = randn(M, 1);

            % Iteration of fields
            for exp = 1:experiments
                A = lambda - sum(J^2./A(neighbours))';
                H = sum(Js.*H(neighbours)./A(neighbours))';
                %H = H/std(H);
                normH(exp) = norm(H);
            end

            % First linear done
            [p, s] = polyfit(experiments-b:experiments, log(normH(experiments-b:end))', 1);
            R2 = 1-s.normr^2/norm(log(normH(experiments-b:end))-mean(log(normH(experiments-b:end))))^2;

            Ps(:, indlambda) = p;
            R2s(indlambda) = R2;
        end

        % Second quadratic fit
        [p, s] = polyfit(Lambda_test, Ps(1, :), 2);
        R2 = 1-s.normr^2/norm(Lambda_test-mean(Lambda_test))^2;
        Rts = roots(p);
        lambdaexp(ens) = Rts(2);%-p(2)/p(1);
    end

    % Record the mean of the obtained eigenvalues for a single value of Delta
    lambda_delta(pt) = mean(lambdaexp);
    disp(pt);
end

clearvars Adj aux count degrees Delta entries experiments exp ind indminus J K kaba l lambda Lambda_test M N n neighbours node p R2 s text Tri
