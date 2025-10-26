
function quiverC3D(x,y,z,u,v,w,varargin)
%quiverC3D creates a 3D quiver plot and adds a color coding. The color coding is
%given by the absolut values of the component vectors. Large values result in colors 
%from the upper end of the used colormap. Plotting parameters have to be changed within 
%the function in this version. Values equal to NaN or inf are set to zero.
%In case of complex arrays the absolute value is used.
% 
%   INPUT:
%       x - 2D matrix, x components of initial points
%       y - 2D matrix, y components of initial points
%       z - 2D matrix, z components of initial points
%       u - 2D matrix, x components of arrows
%       v - 2D matrix, y components of arrows
%       w - 2D matrix, z components of arrows
%       maxNumArrows - a positive integer (non-integer should work as well)
%           number limiting the maximum number of plotted arrows. Since vectors
%           of length zero aren't plotted and the matrices are sampled
%           uniformly, it's possible that fewer arrows will be visible in the
%           end. If maxNumArrows is not given its value is set to inf.
% 
%   OUTPUT:
%       none
% 
%   WARNING!: This function does not create arrows with the same length as
%   would be obtained using quiver3(x, y, z, u, v, w). The implementation
%   requires to disable the automatic scaling. Comparable results are
%   obtained using quiver(x, y, z, u, v, w, 0).
% 
% --------------------------------------------------------------------------------------
% 
%   EXAMPLE:
%       [x,y] = meshgrid(linspace(0,10,100),linspace(0,10,100));
%       z = sin(x.^2 + y.^2);
%       u = exp(-0.2*(x-5).^2 - 0.2*(y-5).^2);
%       v = -u;
%       w = u.*v;
%       quiverC3D(x,y,z,u,v,w,500);
%   
% --------------------------------------------------------------------------------------
% 
%   See also: QUIVER3, LINESPEC, COLORMAP.
% 
%% parse arguments
p = inputParser;
addRequired(p, 'x', @ismatrix);
addRequired(p, 'y', @ismatrix);
addRequired(p, 'z', @ismatrix);
addRequired(p, 'u', @ismatrix);
addRequired(p, 'v', @ismatrix);
addRequired(p, 'w', @ismatrix);
addOptional(p, 'maxNumArrows', inf, @isscalar);
addParameter(p, 'Normalize', false);
addParameter(p, 'LineWidth', 1.5, @isscalar);
addParameter(p, 'maxheadsize', 1, @isscalar);
parse(p, x, y, z, u, v, w, varargin{:});
maxNumArrows = p.Results.maxNumArrows;
normalize = p.Results.Normalize;
lw = p.Results.LineWidth;
hs = p.Results.maxheadsize;
if ~isequal(size(x),size(y),size(z),size(u),size(v),size(w))
    error('X,Y,Z,U,V,W have to be matrices of the same size.');
end
%% initialization
if numel(u) > maxNumArrows
    N = ceil(sqrt(numel(u)/maxNumArrows));
    
    x = x(1:N:end,1:N:end);
    y = y(1:N:end,1:N:end);
    z = z(1:N:end,1:N:end);
    u = u(1:N:end,1:N:end);
    v = v(1:N:end,1:N:end);
    w = w(1:N:end,1:N:end);
end
%% taking care of possible issues
% x = issues(x);
% y = issues(y);
% z = issues(z);
% u = issues(u);
% v = issues(v);
% w = issues(w);
if normalize
    normuvw2 = u.^2 + v.^2 + w.^2;
    normuvw_max = sqrt(max(max(normuvw2)));
    
    u = u./normuvw_max;
    v = v./normuvw_max;
    w = w./normuvw_max;
end
%% colormap definition
C = colormap;
ncolors = size(C, 1);
I = sqrt(u.^2 + v.^2 + w.^2);
% assume that the colormap matrix contains the colors in its rows
Ic = round(I/max(max(I))*ncolors);
Ic(Ic == 0) = 1;
%% plotting
ishold_flag = ishold();
hold on
if numel(u) > ncolors
    % let's take an intelligent approach: group all values which belong to
    % the same color value
    for k = 1:ncolors
        mask = (Ic == k);    
        quiver3(x(mask), y(mask), z(mask), u(mask), v(mask), w(mask), 0, ...
            'Color', C(k, :), 'LineWidth', lw, 'maxheadsize', hs);
    end
else
    % if there are so few values, it will be more efficient to plot each
    % arrow individually
    % linear indexing!
    for k = 1:numel(u)
        quiver3(x(k), y(k), z(k), u(k), v(k), w(k), 0, 'Color', C(Ic(k),:), ...
            'LineWidth', lw, 'maxheadsize', hs);
    end
end
if ~ishold_flag
    hold off
end
end
