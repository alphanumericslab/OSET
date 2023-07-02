classdef HTMLGenerator
    %HTMLGENERATOR A class for generating HTML files
    %
    %   This class generates HTML elements for reporting results
    %
    % Reza Sameni, Sep 2020
    % The open-source electrophysiological toolbox (OSET)
    % www.oset.ir
    %
    properties
        fid % file ID
    end
    
    methods
        function obj = HTMLGenerator(fpath, fname, title, header) % Constructor function
            obj.fid = CreateHTML(fpath, fname, title, header);
        end
        
        function delete(obj) % Destructor function
            CloseHTML(obj.fid);
        end
        
        AddTableHTML(obj, caption);
        AddTableRowHTML(obj, row, style, lastrow);
    end
end

% Creates an HTML file and returns the file id
function fid = CreateHTML(fpath, fname, title, header)
fid = fopen([fpath fname], 'w');
% fprintf(outfid, '<!DOCTYPE html>\n<html lang="en">\n <head>\n\t<meta charset="utf-8">\n\t<title>title</title>\n\t<link rel="stylesheet" href="style.css">\n\t<script src="script.js">\n </head>\n <body>');
fprintf(fid, '<!DOCTYPE html>\n<html lang="en">\n\t<head>\n\t\t<meta charset="utf-8">\n\t\t<title>%s</title>\n', title);
fprintf(fid, '\t\t<style>\n\t\t\ttable { border-collapse: collapse;width: 750%%;}');
fprintf(fid, '\n\t\t\tth {padding: 8px; text-align: left; border-bottom: 1px solid #ddd;}, td {padding: 8px; text-align: left; border-bottom: 1px solid #ddd;}\n\t\t</style>');
fprintf(fid,'\n\t</head>\n\t<body style="font-family: sans-serif">');
fprintf(fid, '\n\t\t<h1>%s</h1>', header);
end

% Close the HTML file
function CloseHTML(fid)
if(fid)
    fprintf(fid, '\n\t</body>\n</html>');
    fclose(fid);
end
end

