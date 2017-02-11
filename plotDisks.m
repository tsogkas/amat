function plotDisks(img,p,r,varargin)
% img: RGB or binary skeleton image
% p,r should be vectors

% Default plotting options
opts.sample = 1;
opts.diskColor = 'r';
opts.diskLineWidth = 0.5;
opts.enhanceVisibility = false;
opts.centerColor = 'g';
opts.centerLineWidth = 1;
opts = parseVarargin(opts, varargin);

if islogical(p)
    p = find(p);
    r = r(p);
end

if opts.sample <= 0 
    error('Sampling factor must be in (0,1]')
elseif opts.sample < 1
    inds = randsample(numel(p), round(opts.sample * numel(p)));
    p = p(inds);
    r = r(inds);
end

[H,W,~] = size(img);
[y,x] = ind2sub([H,W], p);
imshow(img);
% Plot disks
viscircles([x(:),y(:)],r(:),...
            'Color',opts.diskColor,...
            'EnhanceVisibility',opts.enhanceVisibility,...
            'Linewidth',0.5); 

% Plot disk centers
hold on;
plot(x(:),y(:),[opts.centerColor '.'],'Linewidth',opts.centerLineWidth);        
hold off;