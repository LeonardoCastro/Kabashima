
% Random number following a probability distribution p
% The distribution is bounded by kmin and kmax

function r = randomnumber(N, kmin, kmax, kavg, p)
    K = int64([]);

    % Follow the discussion in Appendix A
    % Do an array with floor(Np(k, kavg)) times degree k

    for k = kmin:kmax
        for i = 1:floor(N*p(k, kavg))
            K = [K k];
        end
    end

    % return a random element of the array
    r = K(randi([1, length(K)]));
end
