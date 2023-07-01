function AddTableHTML(obj, caption)
% A member function of HTMLGenerator for starting a table for an HTML table
% Reza Sameni, Sep 2020
% The open-source electrophysiological toolbox (OSET)
% www.oset.ir
%
if(obj.fid)
    fprintf(obj.fid, '\n\t\t<table style="width:750px">\n\t\t\t<caption><h3>%s</h3></caption>', caption);
else
    error('Invalid file ID');
end