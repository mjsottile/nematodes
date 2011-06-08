% seed with a blank row to get started - we clip this off later
lines = [0 0];

% get dimensions of adjacency matrix
[nrows,ncols] = size(adj);

% iterate through adjacency matrix
for i=1:nrows
    % if the current row has a nonzero, do something.  Otherwise, skip.
    if (nnz(adj(i,:))>0)
        % find which columns are nonzero.  this defines the
        % edges in the graph from vertex i to the set of vertices j
        % that it is connected to.
        [~,indices] = find(adj(i,:));
        
        % find the smallest entry in the row, corresponding to the
        % edge of least weight (in this case, weight = euclidean distance)
        [~,minidx] = min(nonzeros(adj(i,:)));
        
        % convert the sparse 1x1 to a regular 1x1 vector
        didx = full(indices(minidx));
        
        % append to the list of valid edges
        lines = [lines; i didx];
    end
end

% clip off the blank row we used to seed the list
lines = lines(2:end,:);


[nlines,~] = size(lines);
figure;
i=1;
plot3([xs(lines(i,1)),xs(lines(i,2))],...
    [ys(lines(i,1)),ys(lines(i,2))],...
    [zs(lines(i,1)),zs(lines(i,2))]);
grid on
axis square
for i=2:nlines
    hold on;
    plot3([xs(lines(i,1)),xs(lines(i,2))],...
        [ys(lines(i,1)),ys(lines(i,2))],...
        [zs(lines(i,1)),zs(lines(i,2))]);
    hold off;
    drawnow;
end
