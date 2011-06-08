lines = [0 0];

[nrows,ncols] = size(adj);
for i=1:nrows
    if (nnz(adj(i,:))>0)
        [~,indices] = find(adj(i,:));
        [~,minidx] = min(nonzeros(adj(i,:)));
        didx = full(indices(minidx));
        lines = [lines; i didx];
    end
end

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
