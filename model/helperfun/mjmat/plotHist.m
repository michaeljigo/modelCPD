% usage:    plotHist(data,binCenters,lineColor,probDist)
% by:       Michael Jigo
% date:     04/27/18
% purpose:  Plot histogram obtained from hist function as a continuous function.
% 
% INPUTS:
% data         vector that will be presented in histogram
% binCenters   a vector defining the center of each bin 
%              (default is determined by hist)
% lineColor    string specifying the color of line plot
% probDist     convert histogram to probability distribution function (default 0)

function figHandle = plotHist(data,binCenters,lineColor,probDist)
if ~exist('lineColor','var') || isempty(lineColor)
   lineColor = 'k';
end

if ~exist('binCenters','var') || isempty(binCenters)
   [~,binCenters] = hist(data);
end

if ~exist('probDist','var') || isempty(probDist)
   probDist = 0;
end

n = sum(~isnan(data(:)));
data = hist(data,binCenters);
if probDist
   data = data./sum(data);
   yLabel = 'probability';
else
   yLabel = '# observations';
end

figHandle = plot(binCenters,data,'Color',lineColor,'Linewidth',4,'linestyle','-');
ylabel(yLabel,'fontname','arial','fontsize',18);
%title(sprintf('n=%i',n));
set(gca,'TickDir','out','box','off','linewidth',2,'fontname','arial','fontsize',12);
