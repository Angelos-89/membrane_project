%-------------------------------------------------------------------------%
TITLENAMES = ["8% pinning","16% pinning","32% pinning"];
FILENAMES = ["pinned_sites_0.txt","pinned_sites_1.txt","pinned_sites_2.txt","pinned_sites_3.txt"];
gridsPerDim = 80; %read this from PARAMS.txt
%--------------------------Data to be plotted-----------------------------%
Zs = {};
for m=1:length(FILENAMES) 
    pinnedSites = importdata(FILENAMES(m));
    dataLen = size(pinnedSites,1);
    xCoor = pinnedSites(:,1);
    yCoor = pinnedSites(:,2);
    xx = 0:1:gridsPerDim-1; yy = xx;
    [X,Y] = meshgrid(xx,yy);
    Z = zeros(gridsPerDim, gridsPerDim);
    for k = 1:dataLen
        i = xCoor(k);
        j = yCoor(k);
        Z(j+1,i+1) = 1; % i,j might be zero
    end 
    Zs{m} = Z;
end
%----------------------------------Now plot-------------------------------%
fig = figure( 'Units', 'normalized', 'Position', [0.0, 0.0, 0.3, 0.5], ...
        'Color', 'white' ) ;
% - Bg axes.
bgAxes = axes( 'Position', [0, 0, 1, 1], 'XColor', 'none', 'YColor', ...
               'none', 'XLim', [0, 1], 'YLim', [0, 1] ) ;
% - Positions.
x = linspace( 0.05, 0.6, 4 ) ;
w = 0.8 * diff( x(1:2) ) ;
y = linspace( 0.07, 0.7, 3 ) ;
h = 0.8 * diff( y(1:2) ) ;

% - Surface 1
axes( 'Position', [0.1, 0.58, 0.35, 0.35] ) ;
s = surf(X,Y,Zs{1}) ;
s.EdgeColor = 'none';
xlim([0,gridsPerDim-1]);
ylim([0,gridsPerDim-1]);
view(2)
ylabel("y (a_0)")
title(TITLENAMES(1))

% - Surface 2
axes( 'Position', [0.55, 0.58, 0.35, 0.35] ) ;
s = surf(X,Y,Zs{2});
s.EdgeColor = 'none';
xlim([0,gridsPerDim-1]);
ylim([0,gridsPerDim-1]);
view(2)
title(TITLENAMES(2))

% - Surface 3
axes( 'Position', [0.1, 0.15, 0.35, 0.35] ) ;
s = surf(X,Y,Zs{3}) ;
s.EdgeColor = 'none';
xlim([0,gridsPerDim-1]);
ylim([0,gridsPerDim-1]);
view(2)
title(TITLENAMES(3))
xlabel("x (a_0)")
ylabel("y (a_0)")

% - Surface 4
axes( 'Position', [0.55, 0.15, 0.35, 0.35] ) ;
s = surf(X,Y,Zs{4});
s.EdgeColor = 'none';
xlim([0,gridsPerDim-1]);
ylim([0,gridsPerDim-1]);
view(2)
title(TITLENAMES(3))
xlabel("x (a_0)")

set( gca, 'Box', 'off' ) ;
saveas(fig,"pinned-distribution.jpg");
