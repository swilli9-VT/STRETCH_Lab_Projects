function SetColorbar(clim, type, cbar_label)
cbar = colorbar;
colormap(jet)
% Dimensions of the colorbar     
% cpos = get(cbar,'position'); 
% cpos(3) = cpos(3)/4 ;   % Reduce the width of colorbar by half
% cpos(2) = cpos(2)+.1 ;
% cpos(4) = cpos(4)-.2 ;
%set(cbar,'Position',cpos) ;
brighten(0.5); 
     
% Title of the colorbar
title(cbar,cbar_label,'fontname','Arial','fontsize',16,'interpreter', 'latex');
%locate = get(cbar,'title');
%tpos = get(locate,'position');
%tpos(3) = tpos(3)+5. ;
%set(locate,'pos',tpos);

% Setting the values on colorbar
%
% get the color limits
numpts = 6;    % Number of points to be displayed on colorbar
if type == "linear"
    caxis(clim);
    ylim(cbar,[clim(1) clim(2)]);
    kssv = linspace(clim(1),clim(2),numpts);
elseif type == "log"
    caxis(10.^clim);
    ylim(cbar,10.^clim);
    kssv = logspace(clim(1),clim(2),numpts);
end
set(cbar,'YtickMode','manual','YTick',kssv); % Set the tickmode to manual
for i = 1:numpts
    imep = num2str(kssv(i),'%+3.2E');
    vasu(i) = {imep} ;
end
set(cbar,'YTickLabel',vasu(1:numpts),'fontsize',10);
%cbar.Layout.Tile = 'south';
cbar.AxisLocation = 'out';
%cbar.Location = "manual";
%cbar.Position = [0.104166666666667,0.11,0.800833333333334,0.026344676180022];

end