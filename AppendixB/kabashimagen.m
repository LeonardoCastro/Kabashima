% Function to generate random grahs following Algorithm 2 of Appendix B.
% The function returns the adjacency matrix

function A = kabashimagen(N, degs)

    % Making sure that the number of degrees is even
    if mod(sum(degs), 2) ~= 0
        auxidx = randi([1, length(degs)]);
        degs(auxidx) = degs(auxidx) + 1;
    end

    % Control in order to assure that the desired distribution of degrees is obtained
    expDegs = zeros(length(degs), 1);
    while ~isequal(expDegs, degs)

        % List of indices given by the degrees
        dsum = sum(degs);
        U = zeros(dsum, 1, 'int64');

        auxidx = 0;
        for i = 1:length(degs)
            for d = 1:degs(i)
                auxidx = auxidx + 1;
                U(auxidx) = i;
            end
        end

        % Allocation of adjacency matrix
        A = spalloc(N, N, dsum);
        clearvars dsum auxidx

        % While there are nodes to link
        while ~isempty(U)

            % If there are not only self-links
            if length( find( U == U(1))) ~= length(U)

                % Choosing two random elements
                i = randi([1, length(U)]);
                j = randi([1, length(U)]);

                % If we are not taking the same element and to avoid self-links
                if i ~= j && U(i) ~= U(j)

                    % If the link already exists, start doing the edge swap
                    if A(U(i), U(j)) == 1 && A(U(j), U(i)) == 1

                        % Find those nodes that are not linked to them
                        [~, J1] = find(~A(U(i), :));
                        [~, J2] = find(~A(U(j), :));

                        [rows, cols] = find( A(J1, J2));

                        if ~isempty(rows)
                            % Delete the existing links

                            link = randi([1, length(rows)]);
                            A(rows(link), cols(link)) = 0;
                            A(cols(link), rows(link)) = 0;

                            % And create the new links
                            A(U(i), rows(link)) = 1;
                            A(rows(link), U(i)) = 1;

                            A(U(j), cols(link)) = 1;
                            A(cols(link), U(j)) = 1;
                        end
                    end


                    % Make link if it does not exist
                    if A(U(i), U(j)) ~= 1 && A(U(j), U(i)) ~=1

                        A(U(i), U(j)) = 1;
                        A(U(j), U(i)) = 1;
                    end


                    % And take the worked nodes out of the list of indices
                    if i < j
                        U(i) = [];
                        U(j-1) = [];
                    end
                    if j < i
                        U(j) = [];
                        U(i-1) = [];
                    end
                end


            % Only self-links left
            elseif length( find( U == U(1))) == length(U)

                 % Find those nodes that are not linked to i
                [~, J1] = find(~A(U(1), :));

                % And choose two that are already linked between them
                [rows, cols] = find( A(J1, J1));

                link1 = randi([1, length(rows)-1]);
                link2 = link1 + 1;

                % Destroy links
                A(rows(link1), cols(link1)) = 0;
                A(cols(link1), rows(link1)) = 0;

                A(rows(link2), cols(link2)) = 0;
                A(cols(link2), rows(link2)) = 0;

                % And create the new links
                A(U(1), rows(link1)) = 1;
                A(rows(link1), U(1)) = 1;

                A(U(2), rows(2)) = 1;
                A(rows(2), U(2)) = 1;

                A(cols(1), cols(2)) = 1;
                A(cols(2), cols(1)) = 1;
                if length(U) > 2
                    U = U(3:end);
                else
                    U = [];
                end
            end
        end
    expDegs = A*ones(N, 1);
    end
end
