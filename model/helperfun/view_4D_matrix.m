% Purpose:  All-purpose viewing function for any 4D matrix. 
%           The last two dimensions will contain the image to be displayed while the first two dimensions correspond to different features/manipulations of the images.

function fig_handle = view_4D_matrix(mat,varargin)

%% Set display parameters
in = {'colormap' 'caxis' 'name' 'x_title' 'y_title' 'x_labels' 'y_labels' 'view_only_y' 'view_only_x'};
val = {'gray' [] '' '' '' '' '' [] []};
disp_params = parseOptionalInputs(in,val,varargin);

% format color axis
if isempty(disp_params.caxis)
   if strcmp(disp_params.colormap,'redblue')
      y_range = [abs(min(mat(:))) max(abs(mat(:)))];
      y_range = [-max(y_range) max(y_range)];
   else
      y_range = [min(mat(:)) max(mat(:))];
   end
   disp_params.caxis = y_range;
end

% format labels
empty_labels = '';
if ~strcmp(disp_params.x_labels,empty_labels)
   if ~ischar(disp_params.x_labels) && ~iscell(disp_params.x_labels)
      disp_params.x_labels = cellfun(@num2str,num2cell(disp_params.x_labels),'uniformoutput',0);
      x_labels = 1;
   end
else
   x_labels = 0;
end

if ~strcmp(disp_params.y_labels,empty_labels)
   if ~ischar(disp_params.y_labels) || ~iscell(disp_params.y_labels)
      disp_params.y_labels = cellfun(@num2str,num2cell(disp_params.y_labels),'uniformoutput',0);
   end
   y_labels = 1;
else
   y_labels = 0;
end

% decide whether to constrain which bands to view
show_y_idx = 1:size(mat,2);
if ~isempty(disp_params.view_only_y)
   show_y_idx = disp_params.view_only_y;
end

show_x_idx = 1:size(mat,1);
if ~isempty(disp_params.view_only_x)
   show_x_idx = disp_params.view_only_x;
end


%% Display matrix
fig_handle = figure('Name',disp_params.name);
plotn = 1;
for y = show_y_idx
   for x = show_x_idx
      subplot(numel(show_y_idx),numel(show_x_idx),plotn);
      imagesc(squeeze(mat(x,y,:,:))); colormap(disp_params.colormap); 
      set(gca,'xtick','','ytick','','box','off','PlotBoxAspectRatio',[1 0.75 0.75]);
      plotn = plotn+1;
     
      if ~strcmp(disp_params.caxis,'none') 
         caxis(disp_params.caxis);
      end
     
     axis square

      % if titles are supplied...
      if y==show_y_idx(end) && x_labels
         % add titles for top-most subplots
         if ~ischar(disp_params.x_labels{x})
            xlabel(sprintf('%.2f',str2num(disp_params.x_labels{x})),'fontname','arial','fontsize',14,'color',[0 0 0]);
         else
            xlabel(disp_params.x_labels{x},'fontname','arial','fontsize',14,'color',[0 0 0]);
         end
      end
      if x==1 && y_labels
         % add titles for left-most subplots
         if ~ischar(disp_params.y_labels{x})
            ylabel(sprintf('%.2f',str2num(disp_params.y_labels{y})),'fontname','arial','fontsize',14,'color',[0 0 0]);
         else
            ylabel(disp_params.y_labels{y},'fontname','arial','fontsize',14,'color',[0 0 0]);
         end
      end 
   end
end
% add titles to horizontal and vertical axes
%mtit('x',disp_params.x_title,'fontsize',16,'color',[0 0 0],'fontweight','bold','fontname','arial');
%mtit('y',disp_params.y_title,'fontsize',16,'color',[0 0 0],'fontweight','bold','fontname','arial');
%set(gcf,'Position',[680 221 555 877],'OuterPosition',[680 221 555 950],'InnerPosition',[680 221 555 877]);
