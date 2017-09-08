function r = randomnumber(N, kmin, kmax, kavg, p)
    K = int64([]);
    
    for k = kmin:kmax
        for i = 1:floor(N*p(k, kavg))
            K = [K k];
        end
    end
    
    r = K(randi([1, length(K)]));
end