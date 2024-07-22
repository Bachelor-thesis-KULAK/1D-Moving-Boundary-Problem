function lgd = add_legend(names,specs1,specs2)
%ADD_LEGEND Add a legend on the current figure
% INPUTS:
%   names: string array of length N, contains the legend names
%   specs1/specs2: either a Nx3 double with RGB of colors, or a string array with linespec strings, or both

ax2 = axes(Position=[0.9,0.92,0,0],Visible='off'); hold(ax2,'on') % axis of size 0
N = length(names); 
lgd_cell = cell(1,N);
if nargin < 3
    if class(specs1) == "string"
        linespec = specs1;
        for i = 1:N, lgd_cell{i} = plot(ax2,nan,linespec(i)); end
    else
        colspec = specs1;
        for i = 1:N, lgd_cell{i} = plot(ax2,nan,Color=colspec(i,:)); end
    end
else
    if class(specs1) == "string"
        linespec = specs1; colspec = specs2;
    else
        colspec = specs1; linespec = specs2;
    end
    for i = 1:N, lgd_cell{i} = plot(ax2,nan,linespec(i),Color=colspec(i,:)); end
end
lgd = legend(ax2,[lgd_cell{:}], names, Interpreter = "latex", fontsize=9);

end
