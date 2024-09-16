% Demonstration of Gaussian convolution and deconvolution.
% Requires Gaussian function 
% Create a rectangular function y, 400 points wide
increment=.01;
%
% Create simulated signal
x=0:increment:20;
y=zeros(size(x));
y(900:1300)=1.3;           

% Add pre-convolution random noise
y=y+.01.*randn(size(y));    % Noise added before the convolution
cw=2; % cw = Convolution width of the physical convolution process (unknown)
widths=[ 1.97 1.98 1.99 2.00 2.01 2.02 2.03];
for trial=1:7,
dw=widths(trial); % dw = Deconvolution width (estimated) must equal cw for perfect results
SmoothWidth=6; % Width of final smoothing to remove high-frequency noise


% Create a centered Gaussian convolution function (maximum at x=zero)
c=gaussian(x,0,cw)+gaussian(x,max(x),cw);  % zero centered Gaussian convolution function
c2=gaussian(x,0,dw)+gaussian(x,max(x),dw);  % zero centered Gaussian deconvolution function

% Convolute rectangular function with Gaussian
yc=ifft(fft(y).*fft(c))./sum(c);  % Create Gaussian convoluted rectangular function

% Add a little bit of post-convolution random noise
yc=yc+.00000001.*randn(size(yc)); % Noise added after the convolution

% Now attempt to recover the original signal by deconvolution (2 methods)
ydc=ifft(fft(yc)./fft(c2)).*sum(c2);  % Deconvolution by fft/ifft
% ydc=deconvgauss(x,yc,w);  % Deconvolution by "deconvgauss" function

% Plot all the steps
subplot(2,2,1); plot(x,y); title('Original signal,y');
subplot(2,2,2); plot(x,yc(1:2001)); title(['yc=y convoluted with Gaussian of width = ' num2str(cw) ]); 
subplot(2,2,3); plot(x,ydc);title(['ydc=Recovered y. Deconvolution width = ' num2str(dw) ]);
subplot(2,2,4); plot(x,fastsmooth(ydc,SmoothWidth,3));title('Smoothed recovered y');
pause(1)
end