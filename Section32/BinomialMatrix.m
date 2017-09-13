%%% Function to obtrain a matrix J where the non-zero entries of the adjacency matrix A
%%% are distributed following p(J) as in Eq. 2.1.

function Adj_modified = BinomialMatrix(Adj, J, Delta)

% Make all positive entries of the Adjacency matrix equal to J
Tri = triu(Adj, 1);
entries = find(Tri);
Tri(entries) = J;

% Then do the entries equal to -J
aux = rand(length(entries), 1);
indminus = find(aux<=(1-Delta)/2);
ind = entries(indminus);

% Return a symmetric matrix
Tri(ind) = -J;
Adj_modified = Tri+Tri';
end
