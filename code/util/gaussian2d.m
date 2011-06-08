function gf = gaussian2d(n)
Dim = [n n];
sig = sqrt(n);
fineness = 1;
[x2d,y2d] = meshgrid(-(Dim(2)-1)/2:fineness:(Dim(2)-1)/2,...
                -(Dim(1)-1)/2:fineness:(Dim(1)-1)/2);
gf    = exp(-(x2d.*x2d + y2d.*y2d)/(2*sig*sig));
gf    = gf/sum(sum(gf))/(fineness^2);
