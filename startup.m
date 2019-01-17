function startup()
% Startup function for AMAT project 

% Create all necessary directories
paths = setPaths();
if ~isdir('./output'), mkdir('output'); end
if ~isdir('./external'), mkdir('external'); end
if ~isdir('./data'), mkdir('data'); end
if ~isdir(paths.amat.models), mkdir(paths.amat.models); end

% Download dependencies
disp('Downloading dependencies...')

% spb-mil -----------------------------------------------------------------
if isdir('./external/spb-mil')  
    disp('-- spb-mil: CHECK!')
else
    disp('Cloning spb-mil...')
    !git clone git@github.com:tsogkas/spb-mil.git external/spb-mil
end
addpath(genpath('external/spb-mil/'))

% matlab-utils ------------------------------------------------------------
if isdir('./external/matlab-utils') 
    disp('-- matlab-utils: CHECK!')
else
    disp('Cloning matlab-utils...')
    !git clone git@github.com:tsogkas/matlab-utils.git external/matlab-utils
end
addpath('external/matlab-utiles')

% L0-smoothing ------------------------------------------------------------
if isdir('./external/L0smoothing')  
    disp('-- L0smoothing: CHECK!')
else
    disp('Downloading L0smoothing code...')
    unzip('http://www.cse.cuhk.edu.hk/leojia/projects/L0smoothing/L0smoothing.zip','external')
end
addpath('external/L0smoothing/code')

% inpaint_nans ------------------------------------------------------------
if isdir('./external/Inpaint_nans')  
    disp('-- L0smoothing: CHECK!')
else
    disp('Downloading L0smoothing code...')
    unzip('http://www.cse.cuhk.edu.hk/leojia/projects/L0smoothing/L0smoothing.zip','external')
end
addpath('external/Inpaint_nans')

% BSDR500 -----------------------------------------------------------------
if isdir('./data/BSR')  
    disp('-- BSDS500: CHECK!')
else
    % untar does not work with this url, so we do it in two steps
    disp('Downloading BSDS500...')
    websave('data/BSR.tgz','http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/BSR/BSR_bsds500.tgz')
    gunzip('data/BSR.tgz'); untar('data/BSR.tgz','data/');
    delete('data/BSR.tgz');
end
disp('-- Added external packages to matlab path: CHECK!')