function [hfig] = plottrajec(posetrajec, fighndl, linemode)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Initializing

n = size(posetrajec,2);
if exist('fighndl','var')
    hfig=fighndl;
else
    hfig = figure();
end

robotsColors = {'-r','-g','-b'};

% Plotting solution estiamtes
hold all
for i=1:n
    pix = posetrajec(13:16:end,i);
    piy = posetrajec(14:16:end,i);
    piz = posetrajec(15:16:end,i);    
    plot3(pix(:).',piy(:).',piz(:).',robotsColors{i},'LineStyle',linemode)    
end
hfig.Children.XLim = [-12,12];
hfig.Children.YLim = [-12,12];
axis square

end

