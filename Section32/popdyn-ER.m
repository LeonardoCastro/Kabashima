%%% Population dynamics algorithm to obtain the distribution of first eigenvalues
%%% for an ERdos-Renyi network

% Initial conditions
c = 4;
J = 1/sqrt(c);
kmin = 0;
kmax = 25;
Delta = 0.;

Np = 1e6;

% Control of sizes
experiments = 40;
b = 25;
ensamble = 500;

ptslambda = 6;

% Allocation of arrays
Lambda_test = linspace(5, 6, ptslambda);
lambdaexp = zeros(ensamble, 1);
Kmax = zeros(ensamble, 1);

for ens = 1:ensamble
    if mod(ens, 100) == 0
        disp(ens);
    end
    meanlnrho = zeros(ptslambda, 1);

    %%%%%%%%% Network %%%%%%%%%
    % Degrees
    degs = randpoissarray(Np+1, kmin, kmax, c);

    kmax_real = max(degs);
    Kmax(ens) = kmax_real;

    % Neighbours
    neighbours = randomneighbours2(degs, Np);

    % Js
    Js = zeros(kmax_real, Np+1);
    ind1 = find(neighbours ~= Np+1);

    Js(ind1) = rand(length(ind1), 1);
    ind2 = find(Js>(1-Delta)/2 );

    Js(ind2) = J;
    Js(Js~=J & Js ~=0) = -J;


    % Non-backtracking extrapolation
    ptlambda = 0;

    % For several values of lambda
    for lambda = Lambda_test
        ptlambda = ptlambda + 1;

        lnrhoexp = zeros(b+1, 1);

        % Allocation of Cavity Fields
        A = rand(Np+1, 1);
        H = rand(Np+1, 1);

        % Allocation of measurements arrays
        normH = zeros(experiments, 1);

        % Iteration of the cavity fields
        for exp = 1:experiments
            A = lambda - sum(J^2./A(neighbours))';
            H = sum(Js.*H(neighbours)./A(neighbours))';
            normH(exp) = norm(H(1:end-1));
        end

        bind = 0;
        for exp = experiments-b:experiments
            bind = bind + 1;
            if normH(exp-1) ~= 0 && normH(exp) ~= 0
                lnrhoexp(bind) = log(normH(exp)/normH(exp-1));
            end
        end

        % Control in case the norm of H goes to zero to quick
        if mean(lnrhoexp) == 0
            for exp = 2:experiments
                if normH(exp-1) ~= 0 && normH(exp) ~= 0
                    lnrhoexp(exp) = log(normH(exp)/normH(exp-1));
                end
            end
        end

        meanlnrho(ptlambda) = mean(lnrhoexp);
    end

    if mean(meanlnrho) ~= 0
        p = polyfit(Lambda_test, meanlnrho', 2);
        Rts = roots(p);
        lambdaexp(ens) = Rts(2);
    end
end

clearvars  A b c Delta degs ens ensamble exp experiments H ind1 ind2 J Js kmax kmax_real kmin lambda neighbours Np pexp ptlambda
