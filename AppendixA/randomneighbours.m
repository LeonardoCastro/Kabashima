function neighbours = randomneighbours(i, deg, Np)

    neighbours = randi([1, Np], deg, length(i));
    control = 0;
    while control == 0 %&& ~isempty(find(neighbours == i, 1)) && ~isempty(find(~diff(sort(neighbours)), 1)) %~isempty(find(~all(diff(neighbours)), 1))
        
        if ~isempty(find(~diff(sort(neighbours)), 1))%~all(diff(neighbours))
            [S, I] = sort(neighbours, 'ascend');
            [row, col] = find(~diff(S));
            neighbours(sub2ind(size(neighbours), I(sub2ind(size(I), row, col)), col)) = randi([1, Np]);
            %ind = find(~diff(sort(neighbours)));
            %neighbours(ind) = randi([1, Np]); %#ok<FNDSB>
            clearvars S I row col
        end
        if ~isempty(find(neighbours == i, 1))
            ind = find(neighbours == i);
            neighbours(ind) = randi([1, Np]); %#ok<FNDSB>
            clearvars ind
        end
        if isempty(find(~diff(sort(neighbours)), 1))
            control = 1;
        end
    end
end