function P=findpeaksT(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% A variant of findpeaksG that measures the peak parameters by constructing
% a triangle around each peak with sides tangent to the sides of the peak.
% Example:
% x=[-10:.2:10];y=exp(-(x+5*randn()).^2)+.1.*randn(size(x));
% P=findpeaksT(x,y,0.003,0.5,7,9,3)
%
% Tom O'Haver, toh@umd.edu, 2016
% Related functions: findpeaksTplot.m, findpeaksG.m
%
P=findpeaksG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype);
sP=size(P);
NumPeaks=sP(1);
% Use the approximate Gaussians fits to locate and fit the peak sides and
% construct triangles upon each peak.
for peak=1:NumPeaks
    center=P(peak,2);
    wid=P(peak,4);
    n0=val2ind(x,center);
    x1=center-wid/2;
    x2=center-wid/6;
    range1=val2ind(x,x1):val2ind(x,x2);
    coeff1=polyfit(x(range1),y(range1),1);
    x4=center+wid/2;
    x3=center+wid/6;
    range2=val2ind(x,x3):val2ind(x,x4);
    coeff2=polyfit(x(range2),y(range2),1);
    x0=(coeff2(2)-coeff1(2))/(coeff1(1)-coeff2(1));
    y0=coeff1(2)+x0*coeff1(1);
    base1=-coeff1(2)/coeff1(1);
    base2=-coeff2(2)/coeff2(1);
    basewidth=base2-base1;
    b=sqrt((center-base1)^2+y0^2);
    a=sqrt((base2-center)^2+y0^2);
    c=basewidth;
    area(peak)=sqrt((a+b-c)*(a-b+c)*(-a+b+c)*(a+b+c))/4;
    P(peak,2)=x0;
    P(peak,3)=y0;
    P(peak,4)=basewidth/2;
    P(peak,5)=area(peak);
end
% ----------------------------------------------------------------------
function P=findpeaksG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% function P=findpeaksG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype)
% Function to locate the positive peaks in a noisy x-y time series data
% set.  Detects peaks by looking for downward zero-crossings
% in the first derivative that exceed SlopeThreshold.
% Returns list (P) containing peak number and position,
% height, width, and area of each peak. Arguments "slopeThreshold",
% "ampThreshold" and "smoothwidth" control peak sensitivity.
% Higher values will neglect smaller features. "Smoothwidth" is
% the width of the smooth applied before peak detection; larger
% values ignore narrow peaks. If smoothwidth=0, no smoothing
% is performed. "Peakgroup" is the number points around the top
% part of the peak that are taken for measurement. If Peakgroup=0
% the local maximum is takes as the peak height and position.
% The argument "smoothtype" determines the smooth algorithm:
%   If smoothtype=1, rectangular (sliding-average or boxcar)
%   If smoothtype=2, triangular (2 passes of sliding-average)
%   If smoothtype=3, pseudo-Gaussian (3 passes of sliding-average)
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% T.C. O'Haver, 1995.  Version 5.2, Last revised September, 2014  
% Fixed bug effecting performance for very small data sets
%
% Examples:
% findpeaksG(0:.01:2,humps(0:.01:2),0,-1,5,5)
% x=[0:.01:50];findpeaks(x,cos(x),0,-1,5,5)
% x=[0:.01:5]';findpeaks(x,x.*sin(x.^2).^2,0,-1,5,5)
% x=[-10:.1:10];y=exp(-(x).^2);findpeaks(x,y,0.005,0.3,3,5,3);
% Find, measure, and plot noisy peak with unknown position
% x=[-10:.2:10];
% y=exp(-(x+5*randn()).^2)+.1.*randn(size(x));
% P=findpeaksG(x,y,0.003,0.5,7,9,3);
% xx=linspace(min(x),max(x));
% yy=P(3).*gaussian(xx,P(2),P(4));
% plot(x,y,'.',xx,yy)
%
% Related functions:
% findvalleys.m, findpeaksL.m, findpeaksb.m, findpeaksb3.m,
% findpeaksplot.m, peakstats.m, findpeaksnr.m, findpeaksGSS.m,
% findpeaksLSS.m, findpeaksfit.m, findsteps.m, findsquarepulse.m, idpeaks.m

% Copyright (c) 2013, 2014 Thomas C. O'Haver
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%
if nargin~=7;smoothtype=1;end  % smoothtype=1 if not specified in argument
if smoothtype>3;smoothtype=3;end
if smoothtype<1;smoothtype=1;end
if smoothwidth<1;smoothwidth=1;end
smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
if smoothwidth>1
    d=fastsmooth(deriv(y),smoothwidth,smoothtype);
else
    d=deriv(y);
end
n=round(peakgroup/2+1);
P=[0 0 0 0 0];
vectorlength=length(y);
peak=1;
AmpTest=AmpThreshold;
for j=2*round(smoothwidth/2)-1:length(y)-smoothwidth-1
    if sign(d(j)) > sign (d(j+1)) % Detects zero-crossing
        if d(j)-d(j+1) > SlopeThreshold % if slope of derivative is larger than SlopeThreshold
            if or(y(j) > AmpTest, y(j+1) > AmpTest)  % if height of peak is larger than AmpThreshold (new version by Anthony Willey)
          % if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold (old version)
                xx=zeros(size(peakgroup));yy=zeros(size(peakgroup));
                for k=1:peakgroup % Create sub-group of points near peak
                    groupindex=j+k-n+2;
                    if groupindex<1, groupindex=1;end
                    if groupindex>vectorlength, groupindex=vectorlength;end
                    xx(k)=x(groupindex);yy(k)=y(groupindex);
                end
                if peakgroup>2
                    [Height, Position, Width]=gaussfit(xx,yy);
                    PeakX=real(Position);   % Compute peak position and height of fitted parabola
                    PeakY=real(Height);
                    MeasuredWidth=real(Width);
                    % if the peak is too narrow for least-squares technique to work
                    % well, just use the max value of y in the sub-group of points near peak.
                else
                    PeakY=max(yy);
                    pindex=val2ind(yy,PeakY);
                    PeakX=xx(pindex(1));
                    MeasuredWidth=0;
                end
                % Construct matrix P. One row for each peak detected,
                % containing the peak number, peak position (x-value) and
                % peak height (y-value). If peak measurement fails and
                % results in NaN, or if the measured peak height is less
                % than AmpThreshold, skip this peak
                if isnan(PeakX) || isnan(PeakY) || PeakY<AmpThreshold
                    % Skip this peak
                else % Otherwiase count this as a valid peak
                    P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth  1.0646.*PeakY*MeasuredWidth];
                    peak=peak+1; % Move on to next peak
                end
            end
        end
    end
end
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
% ----------------------------------------------------------------------
function d=deriv(a)
% First derivative of vector using 2-point central difference.
%  T. C. O'Haver, 1988.
n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1
    d(j)=(a(j+1)-a(j-1)) ./ 2;
end
% ----------------------------------------------------------------------
function SmoothY=fastsmooth(Y,w,type,ends)
% fastbsmooth(Y,w,type,ends) smooths vector Y with smooth
%  of width w. Version 2.0, May 2008.
% The argument "type" determines the smooth type:
%   If type=1, rectangular (sliding-average or boxcar)
%   If type=2, triangular (2 passes of sliding-average)
%   If type=3, pseudo-Gaussian (3 passes of sliding-average)
% The argument "ends" controls how the "ends" of the signal
% (the first w/2 points and the last w/2 points) are handled.
%   If ends=0, the ends are zero.  (In this mode the elapsed
%     time is independent of the smooth width). The fastest.
%   If ends=1, the ends are smoothed with progressively
%     smaller smooths the closer to the end. (In this mode the
%     elapsed time increases with increasing smooth widths).
% fastsmooth(Y,w,type) smooths with ends=0.
% fastsmooth(Y,w) smooths with type=1 and ends=0.
% Example:
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3)= [0 1 4 7 10 7 4 1 1 0]
% fastsmooth([1 1 1 10 10 10 1 1 1 1],3,1,1)= [1 1 4 7 10 7 4 1 1 1]
%  T. C. O'Haver, May, 2008.
if nargin==2, ends=0; type=1; end
if nargin==3, ends=0; end
switch type
    case 1
        SmoothY=sa(Y,w,ends);
    case 2
        SmoothY=sa(sa(Y,w,ends),w,ends);
    case 3
        SmoothY=sa(sa(sa(Y,w,ends),w,ends),w,ends);
end

function SmoothY=sa(Y,smoothwidth,ends)
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
L=length(Y);
for k=1:L-w
    s(k+halfw-1)=SumPoints;
    SumPoints=SumPoints-Y(k);
    SumPoints=SumPoints+Y(k+w);
end
s(k+halfw)=sum(Y(L-w+1:L));
SmoothY=s./w;
% Taper the ends of the signal if ends=1.
if ends==1
    startpoint=(smoothwidth + 1)/2;
    SmoothY(1)=(Y(1)+Y(2))./2;
    for k=2:startpoint
        SmoothY(k)=mean(Y(1:(2*k-1)));
        SmoothY(L-k+1)=mean(Y(L-2*k+2:L));
    end
    SmoothY(L)=(Y(L)+Y(L-1))./2;
end
% ----------------------------------------------------------------------
function [Height, Position, Width]=gaussfit(x,y)
% Converts y-axis to a log scale, fits a parabola
% (quadratic) to the (x,ln(y)) data, then calculates
% the position, width, and height of the
% Gaussian from the three coefficients of the
% quadratic fit.  This is accurate only if the data have
% no baseline offset (that is, trends to zero far off the
% peak) and if there are no zeros or negative values in y.
%
% Example 1: Simplest Gaussian data set
% >> [Height, Position, Width]=gaussfit([1 2 3],[1 2 1]) 
%    returns Height = 2, Position = 2, Width = 2
%
% Example 2: best fit to synthetic noisy Gaussian
% >> x=50:150;y=100.*gaussian(x,100,100)+10.*randn(size(x));
% >> [Height,Position,Width]=gaussfit(x,y) 
%   returns [Height,Position,Width] clustered around 100,100,100.
%
% Example 3: plots data set as circles and best-fit Gaussian as red line
% >> x=[1 2 3 4 5];y=[1 2 2.5 2 1];
% >> [Height,Position,Width]=gaussfit(x,y);
% >> plot(x,y,'o',linspace(0,8),Height.*gaussian(linspace(0,8),Position,Width))
%
% Copyright (c) 2012, Thomas C. O'Haver
maxy=max(y);
for p=1:length(y)
    if y(p)<(maxy/100),y(p)=maxy/100;end
end % for p=1:length(y),
z=log(y);
coef=polyfit(x,z,2);
a=coef(3);
b=coef(2);
c=coef(1);
Height=exp(a-c*(b/(2*c))^2);
Position=-b/(2*c);
Width=2.35482/(sqrt(2)*sqrt(-c));