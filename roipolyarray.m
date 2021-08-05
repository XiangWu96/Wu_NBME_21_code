function d = roipoly(fname)
%ROIPOLY Select polygonal region of interest.
%   Use ROIPOLY to select a polygonal region of interest within an
%   image. ROIPOLY returns a binary image that you can use as a mask for
%   masked filtering.
%
%   BW = ROIPOLY(I,C,R) returns the region of interest selected by the polygon
%   described by vectors C and R. BW is a binary image the same size as I with
%   0's outside the region of interest and 1's inside.
%
%   BW = ROIPOLY(I) displays the image I on the screen and lets you specify
%   the polygon using the mouse. If you omit I, ROIPOLY operates on the image
%   in the current axes. Use normal button clicks to add vertices to the
%   polygon. Pressing <BACKSPACE> or <DELETE> removes the previously selected
%   vertex. A shift-click, right-click, or double-click adds a final vertex to
%   the selection and then starts the fill; pressing <RETURN> finishes the
%   selection without adding a vertex.
%
%   BW = ROIPOLY(x,y,I,xi,yi) uses the vectors x and y to establish a
%   nondefault spatial coordinate system. xi and yi are equal-length vectors
%   that specify polygon vertices as locations in this coordinate system.
%
%   [BW,xi,yi] = ROIPOLY(...) returns the polygon coordinates in xi and
%   yi. Note that ROIPOLY always produces a closed polygon. If the points
%   specified describe a closed polygon (i.e., if the last pair of coordinates
%   is identical to the first pair), the length of xi and yi is equal to the
%   number of points specified. If the points specified do not describe a
%   closed polygon, ROIPOLY adds a final point having the same coordinates as
%   the first point. (In this case the length of xi and yi is one greater than
%   the number of points specified.)
%
%   [x,y,BW,xi,yi] = ROIPOLY(...) returns the XData and YData in x and y; the
%   mask image in BW; and the polygon coordinates in xi and yi.
%
%   If ROIPOLY is called with no output arguments, the resulting image is
%   displayed in a new figure.
%
%   Class Support
%   -------------
%   The input image I can be uint8, uint16, int16, single or double.  The
%   output image BW is logical. All other inputs and outputs are double.
%
%   Remarks
%   -------
%   For any of the ROIPOLY syntaxes, you can replace the input image I with
%   two arguments, M and N, that specify the row and column dimensions of an
%   arbitrary image. If you specify M and N with an interactive form of
%   ROIPOLY, an M-by-N black image is displayed, and you use the mouse to
%   specify a polygon with this image.
%
%   Example
%   -------
%       I = imread('eight.tif');
%       c = [222 272 300 270 221 194];
%       r = [21 21 75 121 121 75];
%       BW = roipoly(I,c,r);
%       figure, imshow(I), figure, imshow(BW)
%
%   See also POLY2MASK, ROIFILT2, ROICOLOR, ROIFILL.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 5.27.4.6 $  $Date: 2004/08/10 01:46:40 $


    % ROIPOLY
    %  Get information from the current figure
    %[x,y,a,hasimage] = getimage;
    xdata=xlim;
    ydata=ylim;
    a=load(fname);
%     if ~hasimage,
%         eid = sprintf('Images:%s:needImageInFigure',mfilename);
%         error(eid,'%s',...
%               'The current figure must contain an image to use ROIPOLY.');
%     end
    [xi,yi] = getline(gcf,'closed'); % Get rect info from the user.
    nrows = size(a,1);
    ncols = size(a,2);


xi = cast_to_double(xi);
yi = cast_to_double(yi);
xdata = cast_to_double(xdata);
ydata = cast_to_double(ydata);
nrows= cast_to_double(nrows);
ncols = cast_to_double(ncols);

%%%
% cast_to_double
%%%
if length(xi)~=length(yi)
    eid = sprintf('Images:%s:xiyiMustBeSameLength',mfilename);
    error(eid,'%s','XI and YI must be the same length.'); 
end

% Make sure polygon is closed.
if (~isempty(xi))
    if ( xi(1) ~= xi(end) || yi(1) ~= yi(end) )
        xi = [xi;xi(1)]; 
        yi = [yi;yi(1)];
    end
end
% Transform xi,yi into pixel coordinates.
roix = axes2pix(ncols, xdata, xi);
roiy = axes2pix(nrows, ydata, yi);

d = poly2mask(roix, roiy, nrows, ncols);



%%%
%%% parse_inputs
%%%




function a = cast_to_double(a)
  if ~isa(a,'double')
    a = double(a);
  end
