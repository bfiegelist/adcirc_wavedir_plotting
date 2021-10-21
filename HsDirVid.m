function HsDirVid(Dirfile_nc,HSfile_nc,meshfile_grd,coord_limits,Hsmax,num_arrows,vidfile_avi,tstepfactor,start_time,my_title)

% Author: Robert Fiegelist (Robert.Fiegelist@uga.edu)

% HSDIRPLOT draws arrows on the current axes for wave direction   
% and plots significant wave height. This is repeated across time and a
% movie file is made. A high-res movie can be made by uncommenting 
% line 110. This MUST be run with acccess to the adcirc_util toolbox made
% by Brian O. Blanton: https://github.com/BrianOBlanton/adcirc_util
% 
%       *Remember to input all file names within quotes
%       INPUT: Dirfile      - .nc file of angles (axis is +x, CCW)
%              HSfile       - .nc file
%              meshfile     - .grd file
%              coord_limits - [longmin longmax latmin latmax]
%              Hsmax        - max Hs value for the colorbar
%              num_arrows   - where num_arrows = n, creates nxn grid over
%                             defined area.
%              vidfile      - Desired name of video file that will be made
%                             in the current directory (must be .avi)
%              tstepfactor  - Time step for video file will be multiplied  
%                             by tstepfactor. tstepfactor = 1 for full
%                             video
%              start_time   - form of "yyyy-MM-DD HR:MIN:SEC" 
%                             Example: '2018-09-23 12:00:00'
%              my_title     - Figure title in quotes
%
%       OUTPUT: .avi movie file
% 
%       NOTES: 1. Grid of arrows may not be returned as perfect square, but 
%              uses nearest locations defined on the mesh.
%             
close all

% Load mesh
fgs = grd_to_opnml(meshfile_grd);

% Load wave direction angles siginficant wave heights, and time
dir_data = ncread(Dirfile_nc,"swan_DIR");
Hs_data =  ncread(HSfile_nc,"swan_HS");
time_data = ncread(HSfile_nc,"time");

% Set up figure and plot colormesh of significant wave height
longmin = coord_limits(1);
longmax = coord_limits(2);
latmin = coord_limits(3);
latmax = coord_limits(4);

% Create grid of desired arrow coordinates and direction 
[x,y] = meshgrid(linspace(longmin,longmax,num_arrows),linspace(latmin,latmax,num_arrows));
x = x(:); y = y(:); arrowdir = zeros(length(x),width(dir_data));

% Create 2-column double of just coordinates
arrow_grdxy = [x.';y.'].';

% Determine number of nodes from the mesh that fall in defined area
count = 0;
for i=1:length(fgs.x)
    if fgs.x(i)<-82 && fgs.x(i)>-90 && fgs.y(i)>25 && fgs.y(i)<31
        count = count+1;
    end
end

% Create and fill new vectors of coordinates only for nodes in mesh that 
% fall within the defined area
fgs2x = zeros(count,1);
fgs2y = zeros(count,1);
fgs2dir = zeros(count,width(dir_data));
j = 1;
for i=1:length(fgs.x)
    if fgs.x(i)<longmax && fgs.x(i)>longmin && fgs.y(i)>latmin && fgs.y(i)<latmax
        fgs2x(j) = fgs.x(i);
        fgs2y(j) = fgs.y(i);
        fgs2dir(j,:) = dir_data(i,:);
        j = j+1;
    end
end
fgs2xy = [fgs2x.';fgs2y.'].';

% Find the nodes that are closest to each point in the nxn grid of arrow 
% coordinates.
closest = zeros(length(arrow_grdxy),2);
smallest = zeros(length(arrow_grdxy),1);

for k = 1:length(arrow_grdxy)
    %compute Euclidean distances for first point:
    distances = sqrt(sum(bsxfun(@minus, fgs2xy,arrow_grdxy(k,:)).^2,2));
    %find the smallest distance and its index:
    smallest(k) = find(distances==min(distances));
    closest(k,:) = fgs2xy(smallest(k),:); % fill coordinates of closest node
end
for m = 1:width(dir_data)
    arrowdir(:,m)=fgs2dir(smallest,m); % find direction angle for that node
end

u = cosd(arrowdir); % cartesian 
v = sind(arrowdir);

% Open video file
vid = VideoWriter(vidfile_avi);
open(vid);

% Set up figure and plot colormesh of significant wave height across
% time.

for m = 1:tstepfactor:width(dir_data)
        figure(1)
        clf
        fh = figure(1);
       % fh.WindowState = 'maximized' % High-res option (large movie file)
        h1 = colormesh2d(fgs,Hs_data(:,m));
        hb = plotbnd(fgs,'Color','k','LineWidth',.5); 
        set(gca, 'XLim', [1 10], 'YLim', [1 10]);
        axis('equal')
        colormap(jet(150))
        colorbar
        caxis([0 Hsmax])
        title(colorbar,'m','FontSize',10);
        xlabel(datestr(time_data(m)/24/3600+datenum(start_time)))
        axis(coord_limits);
        title(my_title)
        
        hold on
     % Plot arrows using closest nodes coordinates and their wave directions
        quiver(closest(:,1),closest(:,2),u(:,m),v(:,m),'Color','black') 
     % Add current figure frame to the video
        frame = getframe(figure(1));
        writeVideo(vid,frame);
end
    close(vid);
end

