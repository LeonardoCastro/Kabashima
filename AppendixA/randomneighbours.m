% Function to obtain an array corresponding to the deg neighbours of the list of node i.

function neighbours = randomneighbours(i, deg, Np)

    % Create array to allocate random neighbours
    neighbours = randi([1, Np], deg, length(i));
    control = 0;

    while control == 0

        % If there are repeated neighbours for a specific node, change it
        if ~isempty(find(~diff(sort(neighbours)), 1))
            [S, I] = sort(neighbours, 'ascend');
            [row, col] = find(~diff(S));
            neighbours(sub2ind(size(neighbours), I(sub2ind(size(I), row, col)), col)) = randi([1, Np]);
            clearvars S I row col
        end

        % If the node is its own neighbour, change it
        if ~isempty(find(neighbours == i, 1))
            ind = find(neighbours == i);
            neighbours(ind) = randi([1, Np]);
            clearvars ind
        end

        % If both tests are passed, then finish
        if isempty(find(~diff(sort(neighbours)), 1))
            control = 1;
        end
    end
end
