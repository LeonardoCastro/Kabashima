% Random array of length N of numbers obtained from a Poissonian distribution
% The Poissonian distribution can be bounded by kmin and kmax

function K = randpoissarray(N, kmin, kmax, kavg)
    % Allocate Array
    K = zeros(N, 1);

    % Compute the number of times that each degree is going to be
    P = floor(N*poisspdf(kmin:kmax, kavg));

    % Control to in order to obtain N degrees
    while sum(P) ~= N

        % If there are less degrees than expected, sum
        if sum(P) < N
           while sum(P) < N
               ind = randi(length(P));
               P(ind) = P(ind) + 1;
           end
        end

        % If there are more degrees than expeced, rest
        if sum(P) > N
           while sum(P) > N
               ind = randi(length(P));
               P(ind) = P(ind) - 1;
           end
        end
    end

    % Do the degree vector
    cum = 0;
    for k = 1:length(P)
        K(cum+1:cum+P(k)) = kmin + k-1;
        cum = cum + P(k);
    end

    % Randomize it
    K = K(randperm(N));
end
