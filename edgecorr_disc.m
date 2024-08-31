function [ec1,ec2] = edgecorr_disc(data,Wx,Wy,varargin)
%%
% This function returns two matrices containing different disc-based 
% edge correction factors for a point pattern in a rectangular observation
% window. The elements of each output array contain the proportion of the
% boundary/area of a disc, centred at a given point and of a particular radius,
% that intersects with the observation window.

% INPUT:
% data      the point pattern of study
% Wx        the x coords of the vertices of the observation window, a rectangle (or square)
% Wy        the y coords of the vertices of the observation window, a rectangle (or square)
% varargin  can contain a fixed distance value, r.

% OUTPUT:
% ec1       a vector or matrix whose entries give the proportion of the 
%           circumference of a disc centred at a given focal point, 
%           potentially with varying radii. 
%           Each row of this array corresponds to a different focal point 
%           and each column to a different radius.
% ec2       a vector or matrix whose entries give the proportion of the
%           area of a disc centred at a given focal point, potentially with 
%           varying radii. Each row of this array corresponds to a distinct
%           focal point 
% 
% last modified by jsmartin.stats@gmail.com in Apr 2016
%%


    n = size(data,1);

% If a fixed distance value r is specified, then we return vectors; if not,
% then each distinct column of the output matrices corresponds to a disc of
% a different radius. In this case, the radii correspond to the distances
% between the focal point and all other points. This is the only difference
% between the two methods below.

    if nargin>3
        r = varargin{1};
        % calculate the distance between the point and the
        % left/right/bottom/top of W
        rx1 = data(:,3)-Wx(1);
        rx2 = Wx(2)-data(:,3);
        ry1 = data(:,4)-Wy(1);
        ry2 = Wy(4)-data(:,4);
        
        % is the distance between the point and the
        % left/right/bottom/top of W greater than r?
        rx1gtdist = rx1>r;
        rx2gtdist = rx2>r;
        ry1gtdist = ry1>r;
        ry2gtdist = ry2>r;
        % is the distance between the point and any of the 4 corners of W
        % less than or equal to r?
        % if the point of interest is within r of any of the corners, then
        % an entire quarter of the boundary will be missing; we can also
        % treat the missing portion of the area of the disc as a separate
        % case
        tooclosecornerx1y1 = rx1.^2 + ry1.^2 <= r^2;
        tooclosecornerx1y2 = rx1.^2 + ry2.^2 <= r^2;
        tooclosecornerx2y1 = rx2.^2 + ry1.^2 <= r^2;
        tooclosecornerx2y2 = rx2.^2 + ry2.^2 <= r^2;
        
        % to calculate the portion that is missing of the area/boundary of
        % the r-disc around the point, we require the angle subtended by
        % the horizontal (vertical) axis and the point on the edge of W
        % below/above (left/right of) this axis that is exactly r from
        % the point of interest.
        % create a vector to store these angles in.
        thetax1 = zeros(n,1);
        thetax2 = zeros(n,1);
        thetay1 = zeros(n,1);
        thetay2 = zeros(n,1);
        % calculate the angle
        thetax1(~rx1gtdist) = acos(rx1(~rx1gtdist)./r);
        thetax2(~rx2gtdist) = acos(rx2(~rx2gtdist)./r);
        thetay1(~ry1gtdist) = acos(ry1(~ry1gtdist)./r);
        thetay2(~ry2gtdist) = acos(ry2(~ry2gtdist)./r);
        
        % boundoverlaparray is used to hold the values calculated for
        % the proportion of the boundary of the disc around each point that
        % intersects the observation window.
        % areaoverlaparray is used to hold the values calculated for the
        % proportion of the area of the disc around each point that
        % intersects the observation window
        % we start by assuming the point is such that its corresponding
        % r-disc is fully contained within W
        boundoverlaparray = 2*pi*r*ones(n,1);
        areaoverlaparray = pi*r*r*ones(n,1);
        % if the point is within r of the left edge but not within r of
        % the bottom left corner of W, subtract the appropriate
        % portion of the circumference of the disc; also, subtract
        % the sector corresponding to the appropriate subtended angle,
        % and then add on the corresponding triangle-section again.
        indmat = (~tooclosecornerx1y1)&(~rx1gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetax1(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetax1(indmat).*(r*r/2) + rx1(indmat).*(r*sin(thetax1(indmat))/2);
        % ...within r of the bottom but not the bottom left corner...
        indmat = (~tooclosecornerx1y1)&(~ry1gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetay1(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetay1(indmat).*(r*r/2) + ry1(indmat).*(r*sin(thetay1(indmat))/2);
        % ...within r of the left but not the top left corner...
        indmat = (~tooclosecornerx1y2)&(~rx1gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetax1(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetax1(indmat).*(r*r/2) + rx1(indmat).*(r*sin(thetax1(indmat))/2);
        % ...within r of the top but not the top left corner...
        indmat = (~tooclosecornerx1y2)&(~ry2gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetay2(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetay2(indmat).*(r*r/2) + ry2(indmat).*(r*sin(thetay2(indmat))/2);
        % ...within r of the right but not the bottom right corner...
        indmat = (~tooclosecornerx2y1)&(~rx2gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetax2(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetax2(indmat).*(r*r/2) + rx2(indmat).*(r*sin(thetax2(indmat))/2);
        % ...within r of the bottom but not the bottom right corner...
        indmat = (~tooclosecornerx2y1)&(~ry1gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetay1(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetay1(indmat).*(r*r/2) + ry1(indmat).*(r*sin(thetay1(indmat))/2);
        % ...within r of the right but not the top right corner...
        indmat = (~tooclosecornerx2y2)&(~rx2gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetax2(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetax2(indmat).*(r*r/2) + rx2(indmat).*(r*sin(thetax2(indmat))/2);
        % ...within r of the top but not the top right corner...
        indmat = (~tooclosecornerx2y2)&(~ry2gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetay2(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetay2(indmat).*(r*r/2) + ry2(indmat).*(r*sin(thetay2(indmat))/2);
        
        % if the point is within r of one of the corners, a quarter of
        % the boundary of the r-disc is removed; similarly, a quarter
        % of the area of the r-disc is removed, and then the
        % corresponding square contained within the intersect of the
        % quarter-circle and W is readded to the area.
        boundoverlaparray(tooclosecornerx1y1) = boundoverlaparray(tooclosecornerx1y1) - pi*r/2;
        areaoverlaparray(tooclosecornerx1y1) = areaoverlaparray(tooclosecornerx1y1) - pi*r*r/4 + rx1(tooclosecornerx1y1).*ry1(tooclosecornerx1y1);
        boundoverlaparray(tooclosecornerx1y2) = boundoverlaparray(tooclosecornerx1y2) - pi*r/2;
        areaoverlaparray(tooclosecornerx1y2) = areaoverlaparray(tooclosecornerx1y2) - pi*r*r/4 + rx1(tooclosecornerx1y2).*ry2(tooclosecornerx1y2);
        boundoverlaparray(tooclosecornerx2y1) = boundoverlaparray(tooclosecornerx2y1) - pi*r/2;
        areaoverlaparray(tooclosecornerx2y1) = areaoverlaparray(tooclosecornerx2y1) - pi*r*r/4 + rx2(tooclosecornerx2y1).*ry1(tooclosecornerx2y1);
        boundoverlaparray(tooclosecornerx2y2) = boundoverlaparray(tooclosecornerx2y2) - pi*r/2;
        areaoverlaparray(tooclosecornerx2y2) = areaoverlaparray(tooclosecornerx2y2) - pi*r*r/4 + rx2(tooclosecornerx2y2).*ry2(tooclosecornerx2y2);

        ec1 = boundoverlaparray./(2*pi*r*ones(n,1));
        ec2 = areaoverlaparray./(pi.*r*r*ones(n,1));
    else
        distsmat = zeros(size(data,1));
        distsmat(tril(true(n),-1)) = pdist(data(:,3:4));
        distsmat = distsmat+distsmat';

        % each row corresponds to a distinct focal point. 
        % calculate the distance between each focal point and the
        % left/right/bottom/top of W:
        rx1 = repmat(data(:,3)-Wx(1),1,n);
        rx2 = repmat(Wx(2)-data(:,3),1,n);
        ry1 = repmat(data(:,4)-Wy(1),1,n);
        ry2 = repmat(Wy(4)-data(:,4),1,n);

        % is the distance between the focal point and the left/right/bottom/top
        % of W greater than the distance between the focal point and any of the
        % other points?
        rx1gtdist = rx1>distsmat;
        rx2gtdist = rx2>distsmat;
        ry1gtdist = ry1>distsmat;
        ry2gtdist = ry2>distsmat;

        % is the distance between the focal point and any of the 4 corners of W
        % less than or equal to the distance between the focal point and any of
        % the other points?
        % if the focal point is within r of any of the corners, then an entire 
        % quarter of the boundary will be missing; we can also treat the 
        % missing portion of the area of the disc as a separate case.
        tooclosecornerx1y1 = rx1.^2 + ry1.^2 <= distsmat.^2;
        tooclosecornerx1y2 = rx1.^2 + ry2.^2 <= distsmat.^2;
        tooclosecornerx2y1 = rx2.^2 + ry1.^2 <= distsmat.^2;
        tooclosecornerx2y2 = rx2.^2 + ry2.^2 <= distsmat.^2;

        % to calculate the portion that is missing of the area/boundary of
        % each disc around the focal point, we require the angle subtended by
        % the horizontal (vertical) axis and the point on the edge of W
        % below/above (left/right of) this axis that is exactly the 
        % corresponding distance from the focal point.
        % create a matrix to store these angles in - each row corresponds to a
        % distinct focal point, each column to the distance from that focal
        % point to another point
        thetax1 = zeros(n,n);
        thetax2 = zeros(n,n);
        thetay1 = zeros(n,n);
        thetay2 = zeros(n,n);
        % calculate the angle
        thetax1(~rx1gtdist) = acos(rx1(~rx1gtdist)./distsmat(~rx1gtdist));
        thetax2(~rx2gtdist) = acos(rx2(~rx2gtdist)./distsmat(~rx2gtdist));
        thetay1(~ry1gtdist) = acos(ry1(~ry1gtdist)./distsmat(~ry1gtdist));
        thetay2(~ry2gtdist) = acos(ry2(~ry2gtdist)./distsmat(~ry2gtdist));

        % boundoverlaparray is used to hold the values calculated for
        % the proportion of the boundaries of the discs of different radii, 
        % centred around each focal point, that intersect the observation window.
        % areaoverlaparray is used to hold the values calculated for the
        % proportion of the area of the discs of different radii, centred 
        % around each focal point that intersect the observation window.
        % We start by assuming the focal points are such that the corresponding
        % discs are fully contained within W
        boundoverlaparray = 2*pi.*distsmat;
        areaoverlaparray = pi.*distsmat.*distsmat;
        % if the focal point is within the relevant distance of the left edge 
        % but not within that distance of the bottom left corner of W, subtract
        % the appropriate portion of the circumference of the disc; also, 
        % subtract the sector corresponding to the appropriate subtended angle,
        % and then add on the corresponding triangle-section again.
        indmat = (~tooclosecornerx1y1)&(~rx1gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetax1(indmat).*distsmat(indmat);
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetax1(indmat).*(distsmat(indmat).*distsmat(indmat)/2) + rx1(indmat).*(distsmat(indmat).*sin(thetax1(indmat))/2);
        % ...within that distance of the bottom but not the bottom left corner...
        indmat = (~tooclosecornerx1y1)&(~ry1gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetay1(indmat).*distsmat(indmat);
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetay1(indmat).*(distsmat(indmat).*distsmat(indmat)/2) + ry1(indmat).*(distsmat(indmat).*sin(thetay1(indmat))/2);
        % ...within that distance of the left but not the top left corner...
        indmat = (~tooclosecornerx1y2)&(~rx1gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetax1(indmat).*distsmat(indmat);
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetax1(indmat).*(distsmat(indmat).*distsmat(indmat)/2) + rx1(indmat).*(distsmat(indmat).*sin(thetax1(indmat))/2);
        % ...within that distance of the top but not the top left corner...
        indmat = (~tooclosecornerx1y2)&(~ry2gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetay2(indmat).*distsmat(indmat);
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetay2(indmat).*(distsmat(indmat).*distsmat(indmat)/2) + ry2(indmat).*(distsmat(indmat).*sin(thetay2(indmat))/2);
        % ...within that distance of the right but not the bottom right corner...
        indmat = (~tooclosecornerx2y1)&(~rx2gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetax2(indmat).*distsmat(indmat);
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetax2(indmat).*(distsmat(indmat).*distsmat(indmat)/2) + rx2(indmat).*(distsmat(indmat).*sin(thetax2(indmat))/2);
        % ...within that distance of the bottom but not the bottom right corner...
        indmat = (~tooclosecornerx2y1)&(~ry1gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetay1(indmat).*distsmat(indmat);
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetay1(indmat).*(distsmat(indmat).*distsmat(indmat)/2) + ry1(indmat).*(distsmat(indmat).*sin(thetay1(indmat))/2);
        % ...within that distance of the right but not the top right corner...
        indmat = (~tooclosecornerx2y2)&(~rx2gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetax2(indmat).*distsmat(indmat);
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetax2(indmat).*(distsmat(indmat).*distsmat(indmat)/2) + rx2(indmat).*(distsmat(indmat).*sin(thetax2(indmat))/2);
        % ...within that distance of the top but not the top right corner...
        indmat = (~tooclosecornerx2y2)&(~ry2gtdist);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - thetay2(indmat).*distsmat(indmat);
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - thetay2(indmat).*(distsmat(indmat).*distsmat(indmat)/2) + ry2(indmat).*(distsmat(indmat).*sin(thetay2(indmat))/2);

        ec1 = boundoverlaparray./(2*pi.*distsmat);
        ec2 = areaoverlaparray./(pi.*distsmat.*distsmat);
    end         

end