function [ec1,ec2] = edgecorr_disc_prlgm(data,Wx,Wy,varargin)
%%
% This function returns two matrices containing different disc-based 
% edge correction factors for a point pattern in a parallelogram observation
% window. The elements of each output array contain the proportion of the
% boundary/area of a disc, centred at a given point and of a particular radius,
% that intersects with the observation window, which is a parallelogram.

% INPUT:
% data      the point pattern of study
% Wx        the x coords of the vertices of the observation window, a parallelogram
% Wy        the y coords of the vertices of the observation window, a parallelogram
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
% last modified by jsmartin.stats@gmail.com in Dec 2017
%%

    n = size(data,1);

% If a fixed distance value r is specified, then we return vectors; if not,
% then each distinct column of the output matrices corresponds to a disc of
% a different radius. In this case, the radii correspond to the distances
% between the focal point and all other points. This is the only difference
% between the two methods below.

    if nargin>2

		centre = [mean(Wx),mean(Wy)];

        % establish the relevant dimensions of the parallelogram by
        % reordering the vertices clockwise about the centre, starting from
        % the first point greater than the centre in both dimensions (i.e.
        % north-east of the centre).
        %...starting from the point with the lowest y-coordinate

        Wleftinds = Wx<centre(1);
        Wrightinds = Wx>=centre(1);
        Wbelowinds = Wy<centre(2);
        Waboveinds = Wy>=centre(2);
        angles = nan(size(Wx));
        angles(Wrightinds&Waboveinds) = atan((Wy(Wrightinds&Waboveinds)-centre(2))./(Wx(Wrightinds&Waboveinds)-centre(1)));
        angles(Wleftinds&Waboveinds) = pi + atan((Wy(Wleftinds&Waboveinds)-centre(2))./(Wx(Wleftinds&Waboveinds)-centre(1)));
        angles(Wleftinds&Wbelowinds) = pi + atan((Wy(Wleftinds&Wbelowinds)-centre(2))./(Wx(Wleftinds&Wbelowinds)-centre(1)));
        angles(Wrightinds&Wbelowinds) = 2*pi + atan((Wy(Wrightinds&Wbelowinds)-centre(2))./(Wx(Wrightinds&Wbelowinds)-centre(1)));
        [~,Worder] = sort(angles);
        firstpoint = find(Wy==min(Wy),1);
        Worder = circshift(Worder,[-find(Worder==firstpoint)+1,0]);
        Wx_cw = Wx(Worder);
        Wy_cw = Wy(Worder);
        sidelengths = sqrt(diff(Wx_cw).^2 + diff(Wy_cw).^2);
        base = sidelengths(1);
        side = sidelengths(2);
        a = min(sidelengths);
        b = max(sidelengths);
        % find the angle between the vectors through the first and second
        % and the second and third points (clockwise)
        AB = [Wx_cw(2)-Wx_cw(1),Wy_cw(2)-Wy_cw(1)];
        BC = [Wx_cw(3)-Wx_cw(2),Wy_cw(3)-Wy_cw(2)];
        skewangle = acos(dot(AB,BC)/(norm(AB)*norm(BC)));
        h = side*sin(skewangle);
        
        %%
        if size(data,2) == 2
            data = [nan(size(data,1),2),data];
        end
        vertex1 = [Wx_cw(1),Wy_cw(1)];
        vertex2 = [Wx_cw(2),Wy_cw(2)];
        vertex3 = [Wx_cw(3),Wy_cw(3)];
        vertex4 = [Wx_cw(4),Wy_cw(4)];
        
        r = varargin{1};
        % calculate the min distance between the point and the edge from
        % window vertices 1to2/2to3/3to4/4to1
        edge12 = repmat(vertex2-vertex1,n,1);
            data_vertex_vec = repmat(vertex2,n,1)-data(:,3:4);
            vec_edge_cross = edge12(:,1).*data_vertex_vec(:,2) - edge12(:,2).*data_vertex_vec(:,1);
                % This is a the absolute value of the cross product. In 3d, 
                % |axb| = |a||b|sin(theta)...since our two vectors a,b lie
                % in the same plane, their cross product will produce a
                % vector in the direction perp to the plane, with abs value
                % equal to the above quantity. Dividing by the abs value of
                % the edge vector, we get the component of the vector from 
                % each data point to vertex2, i.e. the perp distance from
                % the data point to the edge from vertex1 to vertex2
            mindists12 = abs(vec_edge_cross./hypot(edge12(:,1),edge12(:,2)));
        edge23 = repmat(vertex3-vertex2,n,1);
            data_vertex_vec = repmat(vertex3,n,1)-data(:,3:4);
            vec_edge_cross = edge23(:,1).*data_vertex_vec(:,2) - edge23(:,2).*data_vertex_vec(:,1);
            mindists23 = abs(vec_edge_cross./hypot(edge23(:,1),edge23(:,2)));
        edge34 = repmat(vertex4-vertex3,n,1);
            data_vertex_vec = repmat(vertex4,n,1)-data(:,3:4);
            vec_edge_cross = edge34(:,1).*data_vertex_vec(:,2) - edge34(:,2).*data_vertex_vec(:,1);
            mindists34 = abs(vec_edge_cross./hypot(edge34(:,1),edge34(:,2)));
        edge41 = repmat(vertex1-vertex4,n,1);
            data_vertex_vec = repmat(vertex1,n,1)-data(:,3:4);
            vec_edge_cross = edge41(:,1).*data_vertex_vec(:,2) - edge41(:,2).*data_vertex_vec(:,1);
            mindists41 = abs(vec_edge_cross./hypot(edge41(:,1),edge41(:,2)));
        % is the distance between the point and each edge of W greater than
        % r?
        mindists12gtr = mindists12>r;
        mindists23gtr = mindists23>r;
        mindists34gtr = mindists34>r;
        mindists41gtr = mindists41>r;
        % is the distance between the point and any of the 4 corners of W
        % less than or equal to r?
        % if the point of interest is within r of any of the corners, then
        % an entire quarter of the boundary will be missing; we can also
        % treat the missing portion of the area of the disc as a separate
        % case
        mindistsvertex1 = sqrt((data(:,3)-vertex1(1)).^2 + (data(:,4)-vertex1(2)).^2);
        mindistsvertex2 = sqrt((data(:,3)-vertex2(1)).^2 + (data(:,4)-vertex2(2)).^2);
        mindistsvertex3 = sqrt((data(:,3)-vertex3(1)).^2 + (data(:,4)-vertex3(2)).^2);
        mindistsvertex4 = sqrt((data(:,3)-vertex4(1)).^2 + (data(:,4)-vertex4(2)).^2);
        
        tooclosevertex1 =  mindistsvertex1 <= r;
        tooclosevertex2 =  mindistsvertex2 <= r;
        tooclosevertex3 =  mindistsvertex3 <= r;
        tooclosevertex4 =  mindistsvertex4 <= r;
        
        % to calculate the portion that is missing of the area/boundary of
        % the r-disc around the point, we require the angle subtended by
        % the vector normal to the edge that passes through the data point 

        % and the point on the edge of W below/above (left/right of) this 
        % normal vector that is exactly r from the point of interest.
        % Create a vector to store these angles in:
        theta_edge12 = zeros(n,1);
        theta_edge23 = zeros(n,1);
        theta_edge34 = zeros(n,1);
        theta_edge41 = zeros(n,1);
        % calculate the angle
        theta_edge12(~mindists12gtr) = acos(mindists12(~mindists12gtr)./r);
        theta_edge23(~mindists23gtr) = acos(mindists23(~mindists23gtr)./r);
        theta_edge34(~mindists34gtr) = acos(mindists34(~mindists34gtr)./r);
        theta_edge41(~mindists41gtr) = acos(mindists41(~mindists41gtr)./r);

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
        % if the point is within r of the edge between vertices 1 and 2,
        % but not within r of vertex 1, subtract the appropriate portion of
        % the circumference of the disc; also, subtract the sector 
        % corresponding to the appropriate subtended angle, and then add on
        % the corresponding triangle-section again.
        indmat = (~tooclosevertex1)&(~mindists12gtr);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - theta_edge12(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - theta_edge12(indmat).*(r*r/2) + mindists12(indmat).*(r*sin(theta_edge12(indmat)))./2;
        % ...within r of the edge between vertices 1 and 2, but not vertex 2...
        indmat = (~tooclosevertex2)&(~mindists12gtr);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - theta_edge12(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - theta_edge12(indmat).*(r*r/2) + mindists12(indmat).*(r*sin(theta_edge12(indmat)))./2;
        % ...within r of the edge between vertices 2 and 3, but not vertex 2...
        indmat = (~tooclosevertex2)&(~mindists23gtr);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - theta_edge23(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - theta_edge23(indmat).*(r*r/2) + mindists23(indmat).*(r*sin(theta_edge23(indmat)))./2;
        % ...within r of the edge between vertices 2 and 3, but not vertex 3...
        indmat = (~tooclosevertex3)&(~mindists23gtr);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - theta_edge23(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - theta_edge23(indmat).*(r*r/2) + mindists23(indmat).*(r*sin(theta_edge23(indmat)))./2;
        % ...within r of the edge between vertices 3 and 4, but not vertex 3...
        indmat = (~tooclosevertex3)&(~mindists34gtr);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - theta_edge34(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - theta_edge34(indmat).*(r*r/2) + mindists34(indmat).*(r*sin(theta_edge34(indmat)))./2;
        % ...within r of the edge between vertices 3 and 4, but not vertex 4...
        indmat = (~tooclosevertex4)&(~mindists34gtr);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - theta_edge34(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - theta_edge34(indmat).*(r*r/2) + mindists34(indmat).*(r*sin(theta_edge34(indmat)))./2;
        % ...within r of the edge between vertices 4 and 1, but not vertex 4...
        indmat = (~tooclosevertex4)&(~mindists41gtr);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - theta_edge41(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - theta_edge41(indmat).*(r*r/2) + mindists41(indmat).*(r*sin(theta_edge41(indmat)))./2;
        % ...within r of the edge between vertices 4 and 1, but not vertex 1...
        indmat = (~tooclosevertex1)&(~mindists41gtr);
        boundoverlaparray(indmat) = boundoverlaparray(indmat) - theta_edge41(indmat).*r;
        areaoverlaparray(indmat) = areaoverlaparray(indmat) - theta_edge41(indmat).*(r*r/2) + mindists41(indmat).*(r*sin(theta_edge41(indmat)))./2;

        % if the point is within r of one of the vertices, the
        % corresponding arc on the boundary of the r-disc is removed;
        % similarly, the corresponding sector of the r-disc is removed and
        % the required quadrilateral contained within this arc, which still
        % overlaps with W, is re-added
        boundoverlaparray(tooclosevertex1) = boundoverlaparray(tooclosevertex1) - acos(mindists41(tooclosevertex1)./mindistsvertex1(tooclosevertex1)).*r ...
                                                                                - acos(mindists12(tooclosevertex1)./mindistsvertex1(tooclosevertex1)).*r;
        areaoverlaparray(tooclosevertex1) = areaoverlaparray(tooclosevertex1) - acos(mindists41(tooclosevertex1)./mindistsvertex1(tooclosevertex1)).*(r*r/2) + mindists41(tooclosevertex1).*sqrt(mindistsvertex1(tooclosevertex1).^2-mindists41(tooclosevertex1).^2)./2 ...
                                                                              - acos(mindists12(tooclosevertex1)./mindistsvertex1(tooclosevertex1)).*(r*r/2) + mindists12(tooclosevertex1).*sqrt(mindistsvertex1(tooclosevertex1).^2-mindists12(tooclosevertex1).^2)./2;
        boundoverlaparray(tooclosevertex2) = boundoverlaparray(tooclosevertex2) - acos(mindists12(tooclosevertex2)./mindistsvertex2(tooclosevertex2)).*r ...
                                                                                - acos(mindists23(tooclosevertex2)./mindistsvertex2(tooclosevertex2)).*r;
        areaoverlaparray(tooclosevertex2) = areaoverlaparray(tooclosevertex2) - acos(mindists12(tooclosevertex2)./mindistsvertex2(tooclosevertex2)).*(r*r/2) + mindists12(tooclosevertex2).*sqrt(mindistsvertex2(tooclosevertex2).^2-mindists12(tooclosevertex2).^2)./2 ...
                                                                              - acos(mindists23(tooclosevertex2)./mindistsvertex2(tooclosevertex2)).*(r*r/2) + mindists23(tooclosevertex2).*sqrt(mindistsvertex2(tooclosevertex2).^2-mindists23(tooclosevertex2).^2)./2;
        boundoverlaparray(tooclosevertex3) = boundoverlaparray(tooclosevertex3) - acos(mindists23(tooclosevertex3)./mindistsvertex3(tooclosevertex3)).*r ...
                                                                                - acos(mindists34(tooclosevertex3)./mindistsvertex3(tooclosevertex3)).*r;
        areaoverlaparray(tooclosevertex3) = areaoverlaparray(tooclosevertex3) - acos(mindists23(tooclosevertex3)./mindistsvertex3(tooclosevertex3)).*(r*r/2) + mindists23(tooclosevertex3).*sqrt(mindistsvertex3(tooclosevertex3).^2-mindists23(tooclosevertex3).^2)./2 ...
                                                                              - acos(mindists34(tooclosevertex3)./mindistsvertex3(tooclosevertex3)).*(r*r/2) + mindists34(tooclosevertex3).*sqrt(mindistsvertex3(tooclosevertex3).^2-mindists34(tooclosevertex3).^2)./2;
        boundoverlaparray(tooclosevertex4) = boundoverlaparray(tooclosevertex4) - acos(mindists34(tooclosevertex4)./mindistsvertex4(tooclosevertex4)).*r ...
                                                                                - acos(mindists41(tooclosevertex4)./mindistsvertex4(tooclosevertex4)).*r;
        areaoverlaparray(tooclosevertex4) = areaoverlaparray(tooclosevertex4) - acos(mindists34(tooclosevertex4)./mindistsvertex4(tooclosevertex4)).*(r*r/2) + mindists34(tooclosevertex4).*sqrt(mindistsvertex4(tooclosevertex4).^2-mindists34(tooclosevertex4).^2)./2 ...
                                                                              - acos(mindists41(tooclosevertex4)./mindistsvertex4(tooclosevertex4)).*(r*r/2) + mindists41(tooclosevertex4).*sqrt(mindistsvertex4(tooclosevertex4).^2-mindists41(tooclosevertex4).^2)./2;
        
        ec1 = boundoverlaparray./(2*pi*r);
        ec2 = areaoverlaparray./(pi.*r*r);
        
    else
        error('this part of edgecorr_disc_prlgm.m is still under construction');
% % % 
% % % 
% % % 
    end         

end