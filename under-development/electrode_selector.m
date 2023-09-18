function electrode_selector()
    global pointList;  % Declare pointList as global
    
    % Initialize the figure
    hFig = figure('Name', '3D Surface Selector', 'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);
    
    % Initialize pointList
    pointList = [];
    
    % Create axes for 3D plot
    ax3D = axes('Parent', hFig, 'Units', 'normalized', 'Position', [0.05, 0.1, 0.6, 0.8]);
    
    % Create listbox for displaying points
    hListbox = uicontrol('Style', 'listbox', 'Units', 'normalized', 'Position', [0.7, 0.1, 0.25, 0.8]);
    
    % Generate 3D data and plot
    [X, Y] = meshgrid(-5:0.1:5, -5:0.1:5);
    Z = sin(sqrt(X.^2 + Y.^2));
    hSurf = surf(ax3D, X, Y, Z);
    
    % Enable rotation for the 3D plot
    rotate3d on;
    
    % Attach callback for button down on the surface plot
    set(hSurf, 'ButtonDownFcn', @(src, event) recordPoint(src, event, ax3D, hListbox, hFig));
end

function recordPoint(~, ~, ax3D, hListbox, hFig)
    global pointList; % Declare pointList as global
    
    % Check mouse button type
    selectionType = get(hFig, 'SelectionType');
    
    % If right mouse button was clicked
    if strcmp(selectionType, 'alt')
        % Get the clicked point in 3D axes
        clickedPt = get(ax3D, 'CurrentPoint');

        % Extract and store the coordinates
        x = clickedPt(1, 1);
        y = clickedPt(1, 2);
        z = clickedPt(1, 3);

        % Append to pointList
        pointList(end+1, :) = [x, y, z];

        % Update the listbox
        strList = {};
        for i = 1:size(pointList, 1)
            strList{end+1} = ['(', num2str(pointList(i, 1)), ', ', num2str(pointList(i, 2)), ', ', num2str(pointList(i, 3)), ')'];
        end
        set(hListbox, 'String', strList);
    end
end
