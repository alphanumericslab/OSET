function AddTableRowHTML(obj, row, style, lastrow)
% A member function of HTMLGenerator for adding rows to an HTML table
% Reza Sameni, Sep 2020
% The open-source electrophysiological toolbox (OSET)
% www.oset.ir
%
if(obj.fid)
    fprintf(obj.fid, '\n\t\t\t\t<tr>');
    for k = 1 : length(row)
        if(style{k} == 'H')
            fprintf(obj.fid, '<th> %s', row{k});
        elseif(style{k} == 'D')
            fprintf(obj.fid, '<td> %s', row{k});
        end
    end
    if(lastrow)
        fprintf(obj.fid, '\n\t\t</table>');
    end        
else
    error('Invalid file ID');
end