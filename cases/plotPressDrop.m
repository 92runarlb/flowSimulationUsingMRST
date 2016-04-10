function plotPressDrop(G, v1, v2, v3, v4, tit, varargin)

if numel(tit)~= 4
    tit = {'v1', 'v2', 'v3', 'v4'};
end

subplot(2,2,1)
plotCellData(G,v1,varargin{:});
title(tit{1});
axis equal

subplot(2,2,2)
plotCellData(G,v2,varargin{:});
title(tit{2});
axis equal

subplot(2,2,3)
plotCellData(G,v3,varargin{:});
title(tit{3});
axis equal

subplot(2,2,4)
plotCellData(G,v4,varargin{:});
title(tit{4});
axis equal

end
