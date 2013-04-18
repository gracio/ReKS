function [vLineHandleTree hLineHandleTree textHandleTree] = plot2(obj, heightTree, ylbl)
%% PLOT  Plot the tree.

display = 1;    
%% CONSTANTS
    
    LINE_COLOR = [ 0.3 0.3 0.3 ];
    LINE_COLOR_vertical = [ 0.5 0.5 0.5 ];
    
    %% Deal with input
    
    
    if nargin < 2
        heightTree = tree(obj, 1);
    end
    
    if nargin < 3
        ylbl = []; % Will not be plotted
    end

    %% Compute the column width
    
    width = tree(obj, 'clear');

    % Put 1 at the leaves
    iterator = width.depthfirstiterator;
    for i = iterator
       if width.isleaf(i)
           width = width.set(i, 1);
       end
    end
    
    % Cumsum
    width = width.recursivecumfun(@sum);
    
    %% Compute the X *column* position
    
    xcol = tree(width, 'clean');
    
    xcol = xcol.set(1, 1);
    
    previous = 1;
    parent = 1;
    iterator = width.breadthfirstiterator;
    iterator(1) = []; % The root is already done
    
    for i = iterator
       
        newParent = xcol.getparent(i);
        if newParent ~= parent
           % We just changed branch
           parent = newParent;
           previous = xcol.get(parent);
        end
        
        w = width.get(i);
        xcol = xcol.set(i, previous);
        
        previous = previous + w;
    end
    
    %% Compute the actual X position
    
    xpos = tree(width, 'clean');
    iterator = xpos.breadthfirstiterator;
    for i = iterator
        xpos = xpos.set(i, xcol.get(i) + width.get(i)/2 );
    end
    
    % Max of x position
    maxXpos = -1;
    for i = iterator
        xp = xpos.get(i);
        if xp > maxXpos
            maxXpos = xp;
        end
    end
    
    
    
    %% Compute the Y position
    
    ypos = tree(obj, 'clear');
    ypos = ypos.set(1, heightTree.get(1));
    iterator = ypos.depthfirstiterator;
    iterator(1) = []; % Skip the root
    
    maxHeight = heightTree.get(1);
    
    for i = iterator
       parent = ypos.getparent(i);
       parentPos = ypos.get(parent);
       height = heightTree.get(i);
       ypos = ypos.set(i, parentPos + height);
       
       if maxHeight < parentPos + height
           maxHeight = parentPos + height;
       end
    end
    
    %% Prepare the axes
    
    ax = axes(...
        'FontName', 'Courier new', ...
        'FontSize', 9, ...
        'YDir', 'reverse', ...
        'TickDir', 'out', ...
        'XTick', [], ...
        'XTickLabel', '', ...
        'XLim', [0 maxXpos * 1.05], ...
        'XColor', 'w');
    
    if isempty(ylbl)
        set(ax, ...
            'YColor', 'w', ...
            'YTick', [], ...
            'YTickLabel', '')
    else
        ylabel(ylbl, ...
            'HorizontalAlignment', 'right', ...
            'Rotation', 0)
    end
    hold(ax, 'on')
    
    %% A first iteration for the vertical bars
    
    % Prepare holder for the vertical line handles
    vLineHandleTree = tree(obj, 'clear');
    
    iterator = obj.depthfirstiterator;
    for i = iterator
        
        % Vertical bars -> to parent
        
        y1 = ypos.get(i);
        
        if isempty(y1)
            continue
        end
        
        y2 = y1 - heightTree.get(i);
        
        x1 = xpos.get(i);
        x2 = x1;
        
        hl = line([x1 x2], [y1 y2], ...
            'Color', LINE_COLOR_vertical, ...
            'LineWidth', 1);
        
        vLineHandleTree = vLineHandleTree.set(i, hl);
        
    end
    
    % If we were given a height tree, draw white ticks on the tree, a la
    % Tufte.

    if nargin >= 2
        xl = xlim;
        ticks = get(ax, 'YTick');
        
        if nargin >= 3
            % Chop to min and max, a la Tufte
            ticks(end) = maxHeight;
            set(ax, 'YTick', ticks, 'YLim', [0 maxHeight])
        end
        
        % White ticks
        bgColor = get(ax, 'Color');
        for t = ticks
            line( xl + xl(2)/1e2, [t t], ...
                'LineWidth', 2, ...
                'Color', bgColor)
        end
    end
    
    %% New iteration for the bars and the content
        
    % Prepare the holder for the text handles
    textHandleTree = tree(obj, 'clear');
    
    % Prepare the holder for horizontal line handles
    hLineHandleTree = tree(obj, 'clear');
    
    k=1;
    for i = iterator
        
         y1 = ypos.get(i);
        
        if isempty(y1)
            continue
        end
        
        y2 = y1 - heightTree.get(i);
        
        x1 = xpos.get(i);
        
        % The label = content

        content = obj.get(i);
        if ~isempty(content)
            
            if ~ischar(content)
                content = num2str(content);
            end
            
            x12plot(k) = x1;
            y22plot(k) = y2;
            label2plot{k} = {content};
            k=k+1;
            ht = text(x1, y2, {content}, ...  A hack to have text displayed above bars
                'HorizontalAlignment', 'center',...
                'VerticalAlignment', 'middle', ...
                'FontName', 'Courier new', ...
                'Interpreter', 'none', ...
                'FontSize', 12, ...
                'BackgroundColor',[.7 .9 .7]);
            
            textHandleTree = textHandleTree.set(i, ht);
        end

        
        % Horizontal bars -> children
        if obj.isleaf(i)
            continue
        end
        
        children = obj.getchildren(i);
        allX = cell2mat(xpos.Node(children));
        
        y2 = y1;
        x1 = min(allX);
        x2 = max(allX);

        hl = line([x1 x2], [y1 y2], ...
            'Color', LINE_COLOR, ...
            'LineWidth', 5);
        
        hLineHandleTree = hLineHandleTree.set(i, hl);
        
    end
    
    for i=1:length(label2plot)
    ht = text(x12plot(i), y22plot(i), label2plot{i}, ...  A hack to have text displayed above bars
                'HorizontalAlignment', 'center',...
                'VerticalAlignment', 'middle', ...
                'FontName', 'Courier new', ...
                'Interpreter', 'none', ...
                'FontSize', 12, ...
                'BackgroundColor',[.7 .9 .7]);
    end
    
    
end