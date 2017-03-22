function tilefig(FH, nrow, ncol, monitor)
% Tile open figures in a grid on screen.
%   tilefig(), will tile all open figure windows, sorted from top left to
%   bottom right on the bottomleft monitor (if using more than 1 monitor).
%   tilefig(FH) will tile open figure windows with numerical handles FH
%   tilefig(FH,nrow) ( or tilefig([],nrow) ), will plot windows on nrows.
%   tilefig(FH,nrow,ncol), will plot windows in nrow/ncol grid. if grid
%   smaller than the number of handles it will repeat.
%   tilefig(FH,nrow,ncol,monitor), will plot on specified monitor (ordered
%   left to right). if monitor is a vector it will move along the vector
%   using the ncol-nrow grid until all figures are plotted.
%
% INPUT: ---- Feel free to leave inputs empty ([])! ----
% FH      = [vector], figure handles (default: all open figures)
% nrow    = [1x1] Number of rows of figures you want to plot on screen,
%                                      (optional, default Nfigures/Ncol)
% ncol    = [1x1] Number of columns  ,,    ,, (optional, default square)
% monitor = [1x1] 1 or 2; Left (1) or right (2) monitor; (default 1)
%                 rework code if you want top/bottom or 3 monitors.
%          OR [vector]; will move along vector monitors until all figures
%          are plotted.
% Output: ~
% Author: D.A. - DEC.2014 - Raboud University

%% Setting optional inputs
if nargin < 1 || isempty(FH) % No handles -> all open figures
   FH = sort(findobj('Type', 'Figure'));
else
   FH(~ ismember(FH, horzcat(findobj('Type', 'Figure')'))) = [];
end
if nargin < 4 || isempty(monitor)
   monitor = 1;
elseif  any(monitor > size(get(0, 'MonitorPositions'), 1));
   warning('Only one monitor found; Setting "monitor" to 1')
   monitor(monitor > size(get(0, 'MonitorPositions'), 1)) = 1;
end
if nargin < 3 || isempty(ncol)
   if ~ (nargin < 2 || isempty(nrow))
      ncol = ceil(length(FH)/ nrow);
   else
      ncol = ceil(sqrt(length(FH)));
   end
end
if nargin < 2 || isempty(nrow)
   nrow = ceil(length(FH)/ ncol);
end
%% determining figure handles; screensize, figure position ratio
if length(FH) <= 0
   warning('No figures found');
   return
end
sz = sortrows(get(0, 'MonitorPositions'), [1 2]); % monitor 1 = top left
if length(monitor)>1 && length(FH)>nrow*ncol
   functionname = mfilename; % Getting functionname to call it within function
   eval([functionname '(FH(1:(nrow*ncol)),nrow,ncol,monitor(1))'])  
   eval([functionname '(FH((nrow*ncol)+1:end),nrow,ncol,monitor(2:end))']) 
   return;
end
if monitor(1) == 2
   sz(2, 3) = sz(2, 3) - sz(1, 3);
   sz(2, 4) = sz(1, 4);
end
scn_h = min(sz(1:end, 4));
scn_w = min(sz(1:end, 3));
leftside = sz(monitor(1), 1);
if length(FH) > nrow * ncol % if more figures than space --> redo script
   functionname = mfilename; % Getting functionname to call it within function
   eval([functionname '(FH(1:nrow*ncol),nrow,ncol,monitor)'])
   eval([functionname '(FH((nrow*ncol)+1:end),nrow,ncol,monitor)'])
   return
end
fig_h = (scn_h - 50) / nrow;
fig_w = scn_w / ncol;
fig_count = 1; %counter
%% Relocating figures and bringing to front in order;
for i = 1:nrow
   for k = 0:ncol-1
      if fig_count > length(FH)
         return
      end
      set(FH(fig_count), 'OuterPosition', ...  % set position
         [leftside + fig_w * k, scn_h - fig_h * i, fig_w, fig_h]);
      figure(FH(fig_count)); % bring to front;
      fig_count = fig_count + 1; % counter
   end
end

end