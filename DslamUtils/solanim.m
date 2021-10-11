function [] = solanim(ResArray, vidtit, n, BufferSize, perturbedposetrajec, trajecplotmode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Create a video writer object for the output video file and open the object for writing.
% load('starTrajec.mat');
videoname = [vidtit,'.mp4'];

v = VideoWriter(videoname, 'MPEG-4');
v.FrameRate=10;
% ,'FrameRate',4);

open(v);

resnum = size(ResArray,2);
% Generate a set of frames, get the frame from the figure, and then write each frame to the file.
robotsColorWsquare = {'-sr','-sg','-sb'};
robotsOnlyColor = {'r','g','b'};
h1 = figure(1);
for text = 1: resnum %tf 
    
    %Plotting solution up to time t in resnum
    Xhat_text = ResArray{text};
    % Plotting solution estiamtes
    clf(h1)
    hold all   
    pix = Xhat_text.T(13:16:end,:);
    piy = Xhat_text.T(14:16:end,:);
    piz = Xhat_text.T(15:16:end,:);
    for ii=1:n
        plot3(pix(:,ii).',piy(:,ii).',piz(:,ii).',robotsColorWsquare{ii},'MarkerFaceColor',robotsOnlyColor{ii})
    end

    
    % Plotting perturbed and ground truth trajectory
    if exist('perturbedposetrajec','var')
        %%
        for i=1:n
            if trajecplotmode.Noise
                pixPert = perturbedposetrajec(13:16:end,i);
                piyPert = perturbedposetrajec(14:16:end,i);
                pizPert = perturbedposetrajec(15:16:end,i);

                plot3(pixPert(:).',piyPert(:).',pizPert(:).',...
                    robotsOnlyColor{i},'LineWidth',1,'LineStyle','--')
            end
            
            if trajecplotmode.star

                pixStar = starTrajec(13:16:end,i);
                piyStar = starTrajec(14:16:end,i);
                pizStar = starTrajec(15:16:end,i);

                plot3(pixStar(:).',piyStar(:).',pizStar(:).',... 
                    robotsOnlyColor{i},'LineStyle','-','LineWidth',1.5)
            end
            
        end
    end
    axis square
    
    h1.Children.XLim = [-12,12];
    h1.Children.YLim = [-12,12];

   title(sprintf('DSLAM: BufferSize = %i ',BufferSize));
   frame = getframe(h1);
   writeVideo(v,frame);
end

close(v);

end

