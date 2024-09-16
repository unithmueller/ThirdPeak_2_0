function TFitDemo% Keyboard-operated Interactive explorer for measuring absorbance from % instrumentally-broadened transmission profiles. Compares single-% wavelength, Regression, and TFit methods. Simulates photon noise % Random shifts in the 100% T intesity due to background absorption, % and unabsorbed straylight. While the figure window is topmost,% you can adjust these input variables using these pairs of keystrokes:%  True peak A  A/Z   True absorbance of analyte at the peak center, without%                      instrumental broadening, stray light, or noise.%  AbsWidth     S/X  Width of absorption peak (initially = 10)   %  SlitWidth    D/C  Width of instrument function (spectral bandpass) %                      (initially = 10)%  straylight   F/V  Fractional unabsorbed stray light. May be a constant%                     or an array of stray light values at each wavelength. %                     (initially 0.005)%  Noise        G/B  Random noise level (initially 0.05)%  Peak shape     Q  Switch between Gaussian and Lorentzian absorption shape% Other variables:%   ArrayLength = Number of wavelengths (initially 100, set in line 39)  %   y = true transmission spectrum, without noise or broadening%   a = instrument function%   yobsd = noisy instrumentally broadened spectrum%   f = frequency coordinate vector%   x = wavelength coordinate vector%   IzeroShift = Random shifts in the 100% T intesity due to%   backgrouond absorption, cell positioning, dust, etc.% Calculated Outputs:%   Single = Single-wavelengh absorbance at peak center, -log10(T) (line 192)%   WReg. = Absorbance by weighted multiwavelength least-squares regression (line 195)%   TFit = Absorbance by Transmission Fitting (TFit) method (line 201)% Example: just type DemoTFit in the command window.% Tom O'Haver (toh@umd.edu), Version 2.1, November 2011, adds SNR% calculation; W key to Switch between Transmission and Absorbance display. % See http://terpconnect.umd.edu/~toh/spectrum/TFit.htmlglobal x y InstFunction width absorbance SlitWidth ArrayLengthglobal straylight noise NoiseArray IzeroShift shape Tmodeclosefigure(1);clfformat compact% Initial values of the user-adjustable parameters:ArrayLength=100;  % Number of wavelengths (typically 50 - 1000)shape=0; % Absorption lineshape: Shape 0 = Lorentzian; Shape 1 = Gaussian; absorbance=1;  % True absorbance of analyte at peak center (A/Z adjustable)width=ArrayLength/10; % FWHM of absorption peak  (S/X adjustable)SlitWidth=ArrayLength/10; % FWHM of broadening function (D/C adjustable)noise=0.05;  % Random noise level when SlitWidth = 1 (G/B adjustable)straylight=.005; % May be a scalar or a vector of length ArrayLength (F/V adjustable)IzeroShift=noise./20; % Random shifts in the 100% T intesity due to background absorptionTmode=1; % Initial setting of TMode=1 displays transmission spectrum (W to switch)% Calculate and plot graph for the initial conditionsx= [1:ArrayLength]';  % X-axis vector to serve as a wavelength or wavenumber axis.% Perform calculations and plot initial graphs [y,InstFunction]=CalculateTfit(x,width,absorbance,SlitWidth,straylight,noise,NoiseArray,IzeroShift);axis([1 ArrayLength 0 1.1]);disp('Press K key to print list of keystroke commands')% Attaches KeyPress test function to the figure.set(gcf,'KeyPressFcn',@ReadKey)uicontrol('Style','text')% end of outer function% Print list of keyboard commandsdisp('---------------------------------------------------------------')         disp('KEYBOARD COMMANDS')         disp('Peak shape....Q     Toggles between Gaussian and Lorentzian absorption peak shape')         disp('True peak A...A/Z   True absorbance of analyte at peak center, without')         disp('                    instrumental broadening, stray light, or noise.')         disp('AbsWidth......S/X   Width of absorption peak')            disp('SlitWidth.....D/C   Width of instrument function (spectral bandpass)')          disp('Straylight....F/V   Fractional unabsorbed stray light.')          disp('Noise.........G/B   Random noise level')         disp('Re-measure....Spacebar   Re-measure signal')          disp('Switch mode...W     Switch between Transmission and Absorbance display')          disp('Statistics....Tab   Prints table of statistics of 50 repeats')         disp('Cal. Curve....M     Displays analytical calibration curve in Figure 2')         disp('Keys..........K     Print this list of keyboard commands')         function ReadKey(obj,eventdata)global width absorbance SlitWidth global straylight noise shape Tmodekey=get(gcf,'CurrentCharacter');if ischar(key),  switch double(key),     case {122,31}        % When 'z' or up arrow key is pressed, DEcreases "absorbance" by 20%          absorbance=absorbance-.2.*absorbance;          redraw;           case {97,30}        % When 'a' or down arrow key is pressed, INcreases "absorbance" by 20%         absorbance=absorbance+.2.*absorbance;         redraw;     case 120        % When 'x' key is pressed, DEcreases "width" by 1 or 10%          width=width-.1.*width;          redraw;           case 115        % When 's' key is pressed, INcreases "width" by 1 or 10%         width=width+.1.*width;         redraw;     case 99        % When 'c' key is pressed, DEcreases "SlitWidth" by 1 or 10%          SlitWidth=SlitWidth-.1.*SlitWidth;          redraw;           case 100        % When 'd' key is pressed, INcreases "SlitWidth" by 1 or 10%         SlitWidth=SlitWidth+.1.*SlitWidth;         redraw;     case 118        % When 'v' key is pressed, DEcreases "straylight" by 1 or 10%          straylight=straylight-.1.*straylight;          redraw;           case 102        % When 'f' key is pressed, INcreases "straylight" by 1 or 10%         straylight=straylight+.1.*straylight;         redraw;       case 98        % When 'b' key is pressed, DEcreases "noise" by 1 or 10%          noise=noise-.1.*noise;          redraw;           case 103        % When 'g' key is pressed, INcreases "noise" by 1 or 10%         noise=noise+.1.*noise;         redraw;     case 32        % When 'spacebar' key is pressed, re-measures signal         redraw;     case 113        % When 'Q' key is pressed, toggles shape between Gaussian and Lorentzian         if shape==0,             shape=1;         else             shape=0;         end        redraw;     case 119        % When 'W' key is pressed, toggles Tmode  between absorbance and        % tranmission         if Tmode==0,             Tmode=1;         else             Tmode=0;         end        redraw;      case 9        % When 'tab' key is pressed, measures statistics         stats      case 109        % When 'M' key is pressed, computes and displays analaytical curve         absorbancelist=[.01 .02 .03 .05 .1 .2 .3 .5 1 2 3 5 10 20 30 50 100];         CalibrationCurve(absorbancelist,10)      case 107        % When 'K' key is pressed, prints table of keystrokes         disp('KEYBOARD COMMANDS')         disp('Peak shape....Q     Toggles between Gaussian and Lorentzian absorption peak shape')         disp('True peak A...A/Z   True absorbance of analyte at peak center, without')         disp('                    instrumental broadening, stray light, or noise.')         disp('AbsWidth......S/X   Width of absorption peak')            disp('SlitWidth.....D/C   Width of instrument function (spectral bandpass)')          disp('Straylight....F/V   Fractional unabsorbed stray light.')          disp('Noise.........G/B   Random noise level')         disp('Re-measure....Spacebar   Re-measure signal')          disp('Switch mode...W     Switch between Transmission and Absorbance display')          disp('Statistics....Tab   Prints table of statistics of 50 repeats')         disp('Cal. Curve....M     Displays analytical calibration curve in Figure 2')         disp('Keys..........K     Print this list of keyboard commands')               otherwise         UnassignedKey=double(key)       disp('Press k to print out list of keyboard commands')  end % switchend % ischar(key),     function redraw()global x y InstFunction width absorbance SlitWidthglobal straylight noise NoiseArray IzeroShift Tmode% axes(h);ArrayLength=length(x);[y,InstFunction]=CalculateTfit(x,width,absorbance,SlitWidth,straylight,noise,NoiseArray,IzeroShift);function [y,InstFunction]=CalculateTfit(x,width,absorbance,SlitWidth,straylight,noise,NoiseArray,IzeroShift)% Performs calculation and plots graph for DemoTFitglobal z c shape yobsd TrueSpectrum yo Tmode% Define frequency and wavelength coordinate vectorsArrayLength=length(x);j=[-ArrayLength/2:(ArrayLength/2)-1]';f=(ArrayLength/2)-abs(j);% Calculate noisy instrumentally-broadened transmission profile, yobsd,% by convoluting true transmission spectrum y with instrument function InstFunction% and adding noise. if shape,  TrueSpectrum=lorentzian(x,(ArrayLength/2),width);else  TrueSpectrum=gaussian(x,(ArrayLength/2),width);end% Unbroadened transmission profile, plotted as green line on the graph.y=10 .^ (absorbance .* (-TrueSpectrum));fy=fft(y);% Broadening function , plotted as magenta dashed lineInstFunction=gaussian(f,0,SlitWidth);  % defined Gaussian instrument function centered on zerofa=fft(InstFunction);fy1=fy.*fa;             % Convolution by multiplying Fourier transformsyobsd=real(ifft(fy1));  % and inverse transforming the result.yo=yobsd./sum(InstFunction);NoiseArray=randn(size(yo));% ybosd is the broadened transmission profile with noise, plotted as red dots on the graph.% Use only one of the following two lines, to model either photon or% constant noiseObservedNoise=(noise.*NoiseArray).*sqrt(yo).*1/SlitWidth;yobsd=straylight+yo+ObservedNoise;  % Adds simulated photon noise% yobsd=straylight+yo+(noise.*NoiseArray).*1/SlitWidth;  % Adds simulated constant noiseyobsd=yobsd.*(1-straylight);yobsd=yobsd.*(1+IzeroShift.*randn); % Random shifts in Izero figure(1);% Traditional methodsSingleWavelengthAbsorbance=-log10(yobsd(ArrayLength./2));  % Single-wavelength Beer's Law absornanceBackground=ones(size(y));weight=y;WeightedRegression=([weight weight] .* [Background TrueSpectrum])\(-log10(yobsd) .* weight);% TFit methodoptions = optimset('TolX',0.00000001);start=WeightedRegression(2); % Use weighted regression result as the first guessTFitA=fminsearch(@fitM,start,options,yobsd,TrueSpectrum,InstFunction,straylight);if Tmode,    plot(x,real(yobsd),'r.',x,real(y),'g',x,real(c)*z,'b',x,gaussian(x,ArrayLength/2,SlitWidth),'m:');    text(1,1.07,'            Green = Unbroadened T      Magenta = Slit function');    text(1,1.02,'                       Red = Observed T     Blue = Fit to observed T');    text(1,.97,[ '  SNR = '  num2str(SingleWavelengthAbsorbance/std(ObservedNoise))]);    ylabel('Light Intensity/Transmission')    axis([1 ArrayLength 0 1.1]); % Update plotelse    ab=log10(1./real(y));    abf=log10(1./(real(c)*z));    abo=log10(1./real(yobsd));    plot(x,abo,'r.',x,ab,'g',x,abf,'b',x,max(ab).*gaussian(x,ArrayLength/2,SlitWidth),'m:');    text(1,1.07*max(ab),'            Green = Unbroadened A      Magenta = Slit function');     text(1,1.02*max(ab),'            Red = Observed A     Blue = Best fit, converted to A');     text(1,.97*max(ab),[ '  SNR = '  num2str(SingleWavelengthAbsorbance/std(ObservedNoise))]);    ylabel('Absorbance (Press W to switch)')    axis([1 ArrayLength min(abo) 1.1*max(ab)]); % Update plotendif shape,   title([ 'L AbsWidth (S/X) = ' num2str(round(10*width)/10)  '   SlitWidth (D/C)= ' num2str(round(10*SlitWidth)/10) '   Straylight (F/V)= ' num2str(round(1000*mean(straylight))/10) '%   Noise (G/B)= ' num2str(round(1000*noise)/10)  '%' ]);else   title([ 'G AbsWidth (S/X) = ' num2str(round(10*width)/10)  '   SlitWidth (D/C)= ' num2str(round(10*SlitWidth)/10) '   Straylight (F/V)= ' num2str(round(1000*mean(straylight))/10) '%   Noise (G/B)= ' num2str(round(1000*noise)/10)  '%' ]);endabsorbancedisplay=round(10000*absorbance)/10000;SingleWavelengthDisplay=round(10000*real(SingleWavelengthAbsorbance))/10000;Regressiondisplay=round(10000*real(WeightedRegression(2)))/10000;TFitdisplay=round(10000*TFitA(1))/10000;xlabel([ 'True A (A/Z)= ' num2str(absorbancedisplay)  '   Single = ' num2str(SingleWavelengthDisplay) '    WReg. = ' num2str(Regressiondisplay) '  TFit = ' num2str(TFitdisplay)  ])function statsglobal y yo InstFunction absorbance SlitWidth ArrayLengthglobal straylight noise IzeroShift yobsd TrueSpectrumfor k=1:50, % Repeat k times with different random noise samples  yobsd=straylight+yo+((noise/SlitWidth).*randn(size(yo))).*sqrt(yo);   % Add simulated photon noise  yobsd=yobsd.*(1-straylight);  yobsd=yobsd.*(1+IzeroShift.*randn); % Random shifts in Izero from sample to sample after instrument is zeroed  % Conventional methods  SingleWavelengthAbsorbance=-log10(yobsd(ArrayLength./2));  SimpleRegression=TrueSpectrum\(-log10(yobsd));  Background=ones(size(y));  weight=y;  WeightedRegression=([weight weight] .* [Background TrueSpectrum])\(-log10(yobsd) .* weight);  % Curve fitting method  options = optimset('TolX',0.000001);  start=10; % Because of the very large dynamic range of absorbance, two start values are   if SingleWavelengthAbsorbance<1,start=1;end  % used to prevent stalling on local optima.  lam=fminsearch(@fitM,start,options,yobsd,TrueSpectrum,InstFunction,straylight);  results(k,:)=[ SingleWavelengthAbsorbance SimpleRegression WeightedRegression(2) lam];enddisp('------------------------------------------------------------------------')  disp('Statistics of 50 repeat measurements')disp(['True Absorbance = ' num2str(absorbance)])disp('  SingleW    SimpleR    WeightR    TFit') MeanResult=mean(results);% disp(sprintf('%s      %0.4g \t\t%0.4g \t\t%0.4g\t\t%0.4g\t\t%0.4g','Average',MeanResult))% disp(sprintf('%s        %0.2d\t\t%0.2g%% \t\t %0.2g%% \t\t %0.2g%% \t\t %0.2g%%','%RSD',100.*(std(results)./abs(MeanResult))))% disp(sprintf('%s% %0.2d\t\t%0.2g%%\t\t%0.2g%%\t\t%0.2g%%\t\t%0.2g%%','%Accuracy',100.*(MeanResult-absorbance)./absorbance))disp([MeanResult; 100.*(std(results)./abs(MeanResult)); 100.*(MeanResult-absorbance)./absorbance]) disp('The rows are:')disp('Mean value of the measured absorbance')disp('Percent relative standard deviation of the measured absorbance')disp('Percent deviation of the mean from the true absorbance')function CalibrationCurve(absorbancelist,repeats)global x y c z InstFunction width absorbance SlitWidth ArrayLengthglobal straylight noise IzeroShift shape% Compute and plot an analytical calibraiton curve for the absorbances% listed in the array absorbancelist, repeating each measurement "repeats" times.% Define frequency and wavelength coordinate vectorsx=[1:ArrayLength]';j=[-ArrayLength/2:(ArrayLength/2)-1]';f=(ArrayLength/2)-abs(j);% Calculate noisy instrumentally-broadened transmission profile, yobsd,% by convoluting true transmission spectrum y with instrument function% InstFunction and adding noise.  for trial = 1:length(absorbancelist),    absorbance=absorbancelist(trial);    if shape,        TrueSpectrum=lorentzian(x,(ArrayLength/2),width);    else        TrueSpectrum=gaussian(x,(ArrayLength/2),width);    end    y=10 .^ (absorbance .* (-TrueSpectrum));    fy=fft(y);    InstFunction=gaussian(f,0,SlitWidth);  % define Gaussian instrument function centered on zero    fa=fft(InstFunction);    fy1=fy.*fa;            % Convolve the transmission profile with the instrument function (SlitWidth)    yobsd=real(ifft(fy1));  % by multiplying Fourier transforms and inverse transforming the result.    yo=yobsd./sum(InstFunction);    for k=1:repeats, % Repeat k times with different random noise samples        yobsd=straylight+yo+((noise/SlitWidth).*randn(size(yo))).*sqrt(yo);   % Add simulated photon noise        yobsd=yobsd.*(1-straylight);        yobsd=yobsd.*(1+IzeroShift.*randn); % Random shifts in Izero        % Conventional methods        SingleWavelengthAbsorbance=-log10(yobsd(ArrayLength./2));        SimpleRegression=TrueSpectrum\(-log10(yobsd));        Background=ones(size(y));        weight=y;        WeightedRegression=([weight weight] .* [Background TrueSpectrum])\(-log10(yobsd) .* weight);        % Curve fitting method        options = optimset('TolX',0.000001);        start=10; % Because of the very large dynamic range of absorbance, two start values are        if SingleWavelengthAbsorbance<1,start=1;end  % used to prevent stalling on local optima.        lam=fminsearch(@fitM,start,options,yobsd,TrueSpectrum,InstFunction,straylight);        TrueA(trial,k)=[absorbance];SingleWavelength(trial,k)=SingleWavelengthAbsorbance;        SimpleR(trial,k)=SimpleRegression;WeightedR(trial,k)=WeightedRegression(2);        TFit(trial,k)=lam;    end    % (Optional) Plot spectral profiles    figure(2);    plot(x,real(yobsd),'r.',x,real(y),'g',x,real(c)*z,'b',x,gaussian(x,ArrayLength/2,SlitWidth),'m:');    text(5,1.32,'Green = Reference spectrum      Dotted Magenta = Instrument function');    text(5,1.25,'                    Red = Observed T     Blue = Fit to observed T');    xlabel('Wavelength'); ylabel('Transmission');    title([   'True absorbance = ' num2str(absorbance) '    Abs.Width = ' num2str(round(10*width)/10)  '   Inst.Width = ' num2str(round(10*SlitWidth)/10) '    straylight= ' num2str(round(1000*mean(straylight))/10) '%' ]);    axis([0 ArrayLength min(y) 1.1]);    drawnowendloglog(TrueA,TrueA,TrueA,SingleWavelength,'r.',TrueA,TFit,'go',TrueA,SimpleR,'c+',TrueA,WeightedR,'bx')ylim([.01 100]);xlabel('True Peak Absorbance')ylabel('Measured absorbance')title(['    Abs.Width = ' num2str(round(10*width)/10)  '   Inst.Width = ' num2str(round(10*SlitWidth)/10) '    straylight= ' num2str(round(1000*mean(straylight))/10) '%   Noise = ' num2str(round(1000*noise)/10)  '%' ]);text(.02,50,'Red dots = Single wavelength   Cyan + = Simple regression');text(.02,30,' Blue x = Weighted regression   Green o = TFit');function err = fitM(lam,yobsd,Spectra,InstFun,StrayLight)% Fitting function for broadened absorption of any number of components% yobsd =  observed transmission spectrum (column vector)% Sprecta = reference spectra for each component, one component/column% InstFunction = Instrument function or slit function. (column vector)% StrayLight = fractional stray light (scalar or column vector)% Example: % options = optimset('TolX',0.000001);% absorbance=FMINSEARCH('fitM',1,options,[.0123 .0102 .0123 .0147],[.5 1 .5 .2],[1 .5 .0625 .5],.01)% yobsd, Spectra, and InstFunction must have same number of rows (wavelengths)%  T. C. O'Haver, August 2006global zglobal cA = StrayLight + (10 .^ -(Spectra*lam'));fy=fft(A);fa=fft(InstFun);fy1=fy.*fa;                z=real(ifft(fy1))./sum(InstFun);   c = z\yobsd;q = z*c;err = norm(q-yobsd);function g = gaussian(x,pos,wid)%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid%  x may be scalar, vector, or matrix, pos and wid both scalar%  T. C. O'Haver, 1988% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.g = exp(-((x-pos)./(0.6006.*wid)) .^2);function g = lorentzian(x,position,width)% lorentzian(x,position,width) Lorentzian function.% where x may be scalar, vector, or matrix% position and width scalar% T. C. O'Haver, 1988% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);