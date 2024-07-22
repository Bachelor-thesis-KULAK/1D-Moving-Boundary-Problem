function sav(name,kwargs)
%SAV Saves a figure (.fig) and exports to .pdf, in the path "Figures_temp/"
% INPUTS:
%   name: str, name of the file (without extension)
%   optional f: figure handle, the figure to save (default is gcf)
%   optional path: path to filename (default is "Figures_temp/")
%
% OUTPUTS:
%   none
%
% EXAMPLE:
% sav("myfigure")
arguments
    name = "_"
    kwargs.f = gcf;
    kwargs.path = "./Figures/"
end

if ~exist(kwargs.path,"dir"); mkdir(kwargs.path); end

pathtofile = fullfile(kwargs.path, name);
exportgraphics(kwargs.f, pathtofile + ".pdf", ContentType="vector");
savefig(kwargs.f, pathtofile + ".fig", 'compact');

end