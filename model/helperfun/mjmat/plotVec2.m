function [] = plotVec2(M,pairsOfVectors)
% 1) if pairsOfVectors=1, the function will assume that every 2 columns is
% a pair of vectors that will be plotted in the same color. This is useful
% for visualizing the transformations done for question 1c.
% 2) decided to hard-code indices because we've already made sure the input
% will always have 2 rows
% 3) because the the question asked to DRAW the axes, I assumed that I had
% to literally plot lines that would represent the axes and not to restrict
% the figure to have axes that are restricted between -1 and 1

if size(M,1)~=2
    error('Input should have two rows');
end
if ieNotDefined('pairsOfVectors')
    pairsOfVectors = 0; % 1
end

% superficial plotting settings
colors = 'rbgym'; lineStyles = {'-' '--' '-.'};
cCount = 1; lineCount = 1; % initialize color and linestyle counter
legendNames = []; % initialize legend names

% plot 2D vector(s)
figure('Name','vector plot');
for c = 1:size(M,2)
    plotN(c) = plot([0 M(1,c)],[0 M(2,c)],[colors(cCount),lineStyles{lineCount}]); hold on % 2
    scatter(M(1,c),M(2,c),[colors(cCount),'o']);
    % set color of vector
    if pairsOfVectors
        if mod(c,2)==0
            legendNames{end+1} = ['vectorPair ',num2str(c/2)];
            cCount = cCount+1;
        end
    else
        legendNames{end+1} = ['vector ',num2str(c)];
        cCount = cCount+1;
    end
    % reset to initial color if we run out of colors
    if cCount>length(colors)
        cCount = 1;
        lineCount = lineCount+1;
        if lineCount>length(lineStyles)
            warning('Wow you are plotting a lot of vectors...calm down');
            lineCount = 1;
        end
    end
end
plot([-1 1],[0 0],'k-'); plot([0 0],[-1 1],'k-'); % 3
axis equal
% create legend for specific color of vector
if pairsOfVectors
    legend(plotN(2:2:end),legendNames);
else
    legend(plotN,legendNames);
end

end