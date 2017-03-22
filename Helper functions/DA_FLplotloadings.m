function DA_FLplotloadings(LDS,VariableNames,new_color)
% will plot Loadings (LDS varx2) into current plot. column1=x, column2=y.
%
% Input:
%  LDS            = Loadings ( variables x 2 ) where column 1=x axis, column 2=y
%  Variable names = {name1, name2,....}
%  axishandle     = the handle of the axis you want to adjust
% Optional:
%  new_color      = [ N x 3 ] with N=1, all same color or N colors for all
%                  original variables.
% Output
%        none
%
% does not perform 'hold on' or hold off, user should have 'hold on'.
% DA - March 2014

%% set standard parameters
if isa(LDS,'single'); LDS=double(LDS);end % because the <text> can only handle double precision.
axishandle=get(gcf,'CurrentAxes');
if nargin<2 || isempty(VariableNames);
   A=size(LDS,1); VariableNames={['v' num2str(1)]};
   for L1=2:A;
      VariableNames(end+1)={['v' num2str(L1)]};
   end
end
if size(LDS,1)~=length(VariableNames);
   %    disp('ERROR: LDS size - Variable size MISMATCH');
   if size(LDS,1)<length(VariableNames);
      VariableNames=VariableNames(1:size(LDS,1));
   else
      return
   end
end;
arrowcolor=[0.2 .8 0.2];
textcolor=[0.2 .8 0.2] ;
textfontsize=10 ;
textFontWeight='bold';
if nargin>2; arrowcolor=new_color; textcolor=arrowcolor;
else new_color=arrowcolor;
end
%% determining scale
[axisrange] = axis(axishandle);
scaler(1)=(axisrange(1,1))/(-abs(min(LDS(:,1))));
scaler(2)=(axisrange(1,2))/max(LDS(:,1));
scaler(3)=(axisrange(1,3))/(-abs(min(LDS(:,2))));
scaler(4)=(axisrange(1,4))/max(LDS(:,2));
if length(LDS(1,:)) == 3;
    scaler(5)=(axisrange(1,3))/(-abs(min(LDS(:,3))));
    scaler(6)=(axisrange(1,4))/max(LDS(:,3));
end
scaler = min(abs(scaler))*1 ;
LDS_forplot=LDS.*scaler;
%% plottign arrows / text
N_var=size(LDS,1);
if length(LDS(1,:)) <= 2;
    for L1=1:N_var
       if size(new_color,1)>2
          arrowcolor=new_color(L1,:); textcolor=arrowcolor;
       end
       LDS_quiver(L1,:)=quiver(0,0,LDS_forplot(L1,1),LDS_forplot(L1,2),1,'color',arrowcolor, 'linewidth', 2);
       if LDS_forplot(L1,1) > 0
          horzalignment='left';
       else
          horzalignment='right';
       end
       if abs(LDS_forplot(L1,1)/LDS_forplot(L1,2)) > 1
          vertalignment='middle';
       elseif LDS_forplot(L1,2) > 0
          vertalignment='top';
       else
          vertalignment='bottom';
       end
       LDS_text=text(LDS_forplot(L1,1),LDS_forplot(L1,2),VariableNames{L1},'HorizontalAlignment',horzalignment,'VerticalAlignment',vertalignment,'fontsize',textfontsize,'color',textcolor,'FontWeight',textFontWeight) ;
    end
    
elseif length(LDS(1,:)) == 3;
   for L1=1:N_var
       if size(new_color,1)>2
          arrowcolor=new_color(L1,:); textcolor=arrowcolor;
       end
       LDS_quiver(L1,:)=quiver3(0,0,0,LDS_forplot(L1,1),LDS_forplot(L1,2),LDS_forplot(L1,3),1,'color',arrowcolor, 'linewidth', 2);
       if LDS_forplot(L1,1) > 0
          horzalignment='left';
       else
          horzalignment='right';
       end
       if abs(LDS_forplot(L1,1)/LDS_forplot(L1,2)) > 1
          vertalignment='middle';
       elseif LDS_forplot(L1,2) > 0
          vertalignment='top';
       else
          vertalignment='bottom';
       end
       LDS_text=text(LDS_forplot(L1,1),LDS_forplot(L1,2),LDS_forplot(L1,3),VariableNames{L1},'HorizontalAlignment',horzalignment,'VerticalAlignment',vertalignment,'fontsize',textfontsize,'color',textcolor,'FontWeight',textFontWeight) ;
end
end
end