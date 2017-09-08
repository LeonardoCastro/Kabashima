function K = randpoissarray(N, kmin, kmax, kavg)
    K = zeros(N, 1);
    P = floor(N*poisspdf(kmin:kmax, kavg));
    
    while sum(P) ~= N
        if sum(P) < N
           while sum(P) < N
               ind = randi(length(P)); 
               P(ind) = P(ind) + 1;
           end
        end

        if sum(P) > N
           while sum(P) > N
               ind = randi(length(P)); 
               P(ind) = P(ind) - 1;
           end
        end
    end
    cum = 0;
    for k = 1:length(P)
        K(cum+1:cum+P(k)) = kmin + k-1;
        cum = cum + P(k);
    end
    
    K = K(randperm(N));
end