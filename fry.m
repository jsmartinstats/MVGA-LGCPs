function varargout = fry(SPP,varargin)
% This function outputs the Fry plot for the spatial point pattern in SPP
% This is a plot of the vectors that separate all point pairs in SPP
% varargin can contain a minimum and maximum separation distance

% INPUT
% SPP       contains the spatial point pattern of interest
%           This is an n-by-2 matrix containing the x- and y-coords in the first 
%               and second columns, resp.
% varargin  potentially contains:
%               a - the minimum separation distance for contributory point pairs
%               b - the maximum separation distance for contributory point pairs
% 
% OUTPUT
% varargout contains an m-by-2 matrix with the x- and y-coords for the Fry points
%
% last modified by jsmartin.stats@gmail.com in Nov 2016

    if nargin>2
        a = varargin{1};
        b = varargin{2};
    elseif nargin>1
        a = 0;
        b = varargin{1};
    else
        a = 0;
        b = realmax;
    end
        

    n = size(SPP,1);

    % % use basic trig to calculate the angle between two points - this is 
    % % the angle between the line from A to B and the x-axis
    point1 = nan(n,n,2);
    point2 = nan(n,n,2);
    point1(:,:,1) = repmat(SPP(:,1),1,n);
    point1(:,:,2) = repmat(SPP(:,2),1,n);
    point2(:,:,1) = repmat(SPP(:,1)',n,1);
    point2(:,:,2) = repmat(SPP(:,2)',n,1);

    % % below names are according to their interpretation when the connecting 
    % % vector is in the first quadrant
    opp = point2(:,:,2)-point1(:,:,2);
    adj = point2(:,:,1)-point1(:,:,1);

    pdists = squareform(pdist(SPP));
    ysepdist = opp(~logical(eye(n)) & pdists>=a & pdists<=b); % we do not consider the diagonal elements of pdists; this corresponds
    xsepdist = adj(~logical(eye(n)) & pdists>=a & pdists<=b); % to each point's distance to itself...and is irrelevant!

    varargout{1} = [xsepdist,ysepdist];
end