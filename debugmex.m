img = imresize(imread('google.jpg'), 0.5);
disp('Loaded image')
mat = AMAT();
mat.initialize(img);
mat.computeEncodings();
mat.computeCosts();
disp('Lets see if this thing is working!')
% dbstop in AMAT at 589; % line of entry for mexFile
dbmex on;
mat.setCoverMex();