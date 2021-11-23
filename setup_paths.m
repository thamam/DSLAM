% Add this directory and all subdirectories to path for session.
addpath(genpath(fileparts(which(mfilename))));

% Remove assets from path.
rmpath(genpath('./assets/'));
rmpath(genpath('./Version_4/'));

importmanopt;