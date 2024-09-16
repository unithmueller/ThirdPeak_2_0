function [exportParamNameArray, exportParamArray, nBins, minv, maxv] = plotJDDHistogramAndEverythingElse(Axes, jumpDistances, Nb, dimension, param, beta, betaBoot, modelProb,strucNames, nameArray, indArray,diffModes,tau, dr, ri, isPixel, lengthUnit, pxsize, timeunit, analyze3dAs2d, Ntracks)
%Function to plot the results from the anaylsis of JDDs

%% Determine dimensionality
if size(dimension,2) == 1
    dimVal = 1;
elseif size(dimension,2) == 2
    dimVal = 2;
elseif size(dimension,2) == 3
    dimVal = 3;
end
if analyze3dAs2d
    dimVal = 2;
end

%% Plot the histogram
jddHist = histogram(Axes, jumpDistances, Nb);
binLimits = jddHist.BinLimits;
xlim(Axes, [binLimits(1), binLimits(2)]);
nBins = Nb;
%currentBinEdges = jddHist.BinEdges;
minv = binLimits(1);
maxv = binLimits(2);
hold(Axes, "on");

if isPixel
    xlabel(Axes, 'Jump Distance [px]');
else
    xlabel(Axes, sprintf('Jump Distance [%s]',lengthUnit));
end
ylabel(Axes, 'Count')
if size(dimension,2) == 1
    if dimension == 3
        title(Axes, 'Jump Distance Distribution in 1D for X');
    elseif dimension == 4
        title(Axes, 'Jump Distance Distribution in 1D for Y');
    elseif dimension == 5
        title(Axes, 'Jump Distance Distribution in 1D for Z');
    end
    dimVal = 1;
elseif size(dimension,2) == 2
    title(Axes, 'Jump Distance Distribution in 2D')
elseif size(dimension,2) == 3
    if analyze3dAs2d
        title(Axes, 'Jump Distance Distribution of 3D data in 2D')
    else
        title(Axes, 'Jump Distance Distribution in 3D')
    end
end

%% Plotting parameters
linWidth = 1.5;
colors = hsv(14);

%Divide first by dimensionality, then by diffusion type
%1 = FreeDiff   6 = Free+Directed      11 = Directed+Confined
%2 = Directed   7 = Free+Anomal        12 = Anomal+Anomal
%3 = Anomalous  8 = Free+Confined      13 = Anomal+Confined
%4 = Confined   9 = Directed+Directed  14 = Confined+Confined
%5 = Free+Free  10 = Directed+Anomal
legendLabelArray = strings(size(diffModes,1),1);
exportParamNameArray = cell(size(diffModes,1),2);
exportParamArray = cell(size(diffModes,1),3);
%% generate the y data
for i = 1:size(diffModes,1)

    switch diffModes(i)
        case 1
            yval = ND(beta(1,1),dr,ri,tau, dimVal, Ntracks);
            col = colors(1,:);
            paramErrors = betaBoot(1);
            locModelProb = modelProb(modelProb(:,1) == 1,2);
            exportParamNameArray{i,1} = 1;
            exportParamNameArray{i,2} = "Free Diffusion D";
            exportParamArray{i,1} = beta(1,1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Free Diffusion, D=',num2str(beta(1,1)), "±", num2str(betaBoot(1,1)),"; Prob:",num2str(locModelProb)]);
        case 2
            yval = NV(beta(2:3,1), dr,ri, tau, dimVal, Ntracks);
            col = colors(2,:);
            paramErrors = betaBoot(2:3,1);
            locModelProb = modelProb(modelProb(:,1) == 2,2);
            exportParamNameArray{i,1} = 2;
            exportParamNameArray{i,2} = ["Directed Motion V"; "Directed Motion D"];
            exportParamArray{i,1} = beta(2:3,1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Directed Motion, V=',num2str(beta(2,1)), "±", num2str(betaBoot(2,1))," D=,",num2str(beta(3,1)),"±", num2str(betaBoot(3,1)), "; Prob:",num2str(locModelProb)]);
        case 3
            yval = NA(beta(4:5,1), dr,ri, tau, dimVal, Ntracks);
            col = colors(3,:);
            paramErrors = betaBoot(4:5,1);
            locModelProb = modelProb(modelProb(:,1) == 3,2);
            exportParamNameArray{i,1} = 3;
            exportParamNameArray{i,2} = ["Anomal Motion D"; "Anomal Motion alpha"];
            exportParamArray{i,1} = beta(4:5,1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Anomal Motion, D=',num2str(beta(4,1)), "±", num2str(betaBoot(4,1))," A=,",num2str(beta(5,1)),"±", num2str(betaBoot(5,1)), "; Prob:",num2str(locModelProb)]);
        case 4
            yval = NConf(beta(6:7,1), dr,ri, tau, dimVal, Ntracks);
            col = colors(4,:);
            paramErrors = betaBoot(6:7,1);
            locModelProb = modelProb(modelProb(:,1) == 4,2);
            exportParamNameArray{i,1} = 4;
            exportParamNameArray{i,2} = "Confined Motion D";
            exportParamArray{i,1} = beta(6:7,1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Confined Motion, D=',num2str(beta(6,1)), "±", num2str(betaBoot(6,1)), "; Prob:",num2str(locModelProb)]);
        case 5
            yval = NDD(beta(8:10,1), dr,ri, tau, dimVal, Ntracks);
            col = colors(5,:);
            paramErrors = betaBoot(8:10,1);
            locModelProb = modelProb(modelProb(:,1) == 5,2);
            exportParamNameArray{i,1} = 5;
            exportParamNameArray{i,2} = ["Free+Free D1";"Free+Free D2";"Free+Free FD"];
            exportParamArray{i,1} = beta(8:10,1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Free+Free Diffusion, D1=',num2str(beta(8,1)), "±", num2str(betaBoot(8,1))," D2=,",num2str(beta(9,1)),"±", num2str(betaBoot(9,1))," FD=,",num2str(beta(10,1)),"±", num2str(betaBoot(10,1)), "; Prob:",num2str(locModelProb)]);
        case 6
            yval = NDV(beta(11:14,1), dr,ri, tau, dimVal, Ntracks);
            col = colors(6,:);
            paramErrors = betaBoot(11:14,1);
            locModelProb = modelProb(modelProb(:,1) == 6,2);
            exportParamNameArray{i,1} = 6;
            exportParamNameArray{i,2} = ["Free+Directed Df"; "Free+Directed V";"Free+Directed Dv";"Free+Directed FD"];
            exportParamArray{i,1} = beta(11:14,1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Free+Directed Diffusion, D=',num2str(beta(11,1)), "±", num2str(betaBoot(11,1))," V=,",num2str(beta(12,1)),"±", num2str(betaBoot(12,1))," Dv=,",num2str(beta(13,1)),"±", num2str(betaBoot(13,1))," FD=,",num2str(beta(14,1)),"±", num2str(betaBoot(14,1)), "; Prob:",num2str(locModelProb)]);
        case 7
            yval = NDA(beta(15:18,1), dr,ri, tau, dimVal, Ntracks);
            col = colors(7,:);
            paramErrors = betaBoot(15:18,1);
            locModelProb = modelProb(modelProb(:,1) == 7,2);
            exportParamNameArray{i,1} = 7;
            exportParamNameArray{i,2} = ["Free+Anomal Df"; "Free+Anomal Da"; "Free+Anomal alpha"; "Free+Anomal FD"];
            exportParamArray{i,1} = beta(15:18,1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Free+Anomal Diffusion, D1=',num2str(beta(15,1)), "±", num2str(betaBoot(15,1))," Da=,",num2str(beta(16,1)),"±", num2str(betaBoot(16,1))," A=,",num2str(beta(17,1)),"±", num2str(betaBoot(17,1))," FD=,",num2str(beta(18,1)),"±", num2str(betaBoot(18,1)), "; Prob:",num2str(locModelProb)]);
        case 8
            yval = NDConf(beta([19,20,22],1), dr,ri, tau, dimVal, Ntracks);
            col = colors(8,:);
            paramErrors = betaBoot([19,20,22],1);
            locModelProb = modelProb(modelProb(:,1) == 8,2);
            exportParamNameArray{i,1} = 8;
            exportParamNameArray{i,2} = ["Free+Confined Df"; "Free+Confined Dc"; "Free+Confined FD"];
            exportParamArray{i,1} = beta([19,20,22],1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Free+Confined Diffusion, D1=',num2str(beta(19,1)), "±", num2str(betaBoot(19,1))," Dc=,",num2str(beta(20,1)),"±", num2str(betaBoot(20,1))," FD=,",num2str(beta(22,1)),"±", num2str(betaBoot(22,1)), "; Prob:",num2str(locModelProb)]);
        case 9
            yval = NVV(beta([25,23,26,24,27],1), dr,ri, tau, dimVal, Ntracks);
            col = colors(9,:);
            paramErrors = betaBoot([25,23,26,24,27],1);
            locModelProb = modelProb(modelProb(:,1) == 9,2);
            exportParamNameArray{i,1} = 9;
            exportParamNameArray{i,2} = ["Directed+Directed V1";"Directed+Directed D1";"Directed+Directed V2";"Directed+Directed D2";"Directed+Directed FD"];
            exportParamArray{i,1} = beta([25,23,26,24,27],1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Directed+Directed Motion, V1=',num2str(beta(25,1)), "±", num2str(betaBoot(25,1))," D1=,",num2str(beta(23,1)),"±", num2str(betaBoot(23,1)),' V2=',num2str(beta(26,1)), "±", num2str(betaBoot(26,1))," D2=,",num2str(beta(24,1)),"±", num2str(betaBoot(24,1))," FD=,",num2str(beta(27,1)),"±", num2str(betaBoot(27,1)), "; Prob:",num2str(locModelProb)]);
        case 10
            yval = NDirAno(beta([29,28,30,31,32],1), dr,ri, tau, dimVal, Ntracks);
            col = colors(10,:);
            paramErrors = betaBoot([29,28,30,31,32],1);
            locModelProb = modelProb(modelProb(:,1) == 10,2);
            exportParamNameArray{i,1} = 10;
            exportParamNameArray{i,2} = ["Directed+Anomal V";"Directed+Anomal Dv";"Directed+Anomal Da";"Directed+Anomal alpha";"Directed+Anomal FD"];
            exportParamArray{i,1} = beta([29,28,30,31,32],1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Directed+Anomal Diffusion, V1=',num2str(beta(29,1)), "±", num2str(betaBoot(29,1))," Dv=,",num2str(beta(28,1)),"±", num2str(betaBoot(28,1))," Da=,",num2str(beta(30,1)),"±", num2str(betaBoot(30,1))," A=,",num2str(beta(31,1)),"±", num2str(betaBoot(31,1))," FD=,",num2str(beta(32,1)),"±", num2str(betaBoot(32,1)), "; Prob:",num2str(locModelProb)]);
        case 11
            yval = NDirConf(beta([34,33,35,37],1), dr,ri, tau, dimVal, Ntracks);
            col = colors(11,:);
            paramErrors = betaBoot([34,33,35,37],1);
            locModelProb = modelProb(modelProb(:,1) == 11,2);
            exportParamNameArray{i,1} = 11;
            exportParamNameArray{i,2} = ["Directed+Confined V";"Directed+Confined Dv";"Directed+Confined Dc";"Directed+Confined FD"];
            exportParamArray{i,1} = beta([34,33,35,37],1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Directed+Confined Motion, V1=',num2str(beta(34,1)), "±", num2str(betaBoot(34,1))," Dv=,",num2str(beta(33,1)),"±", num2str(betaBoot(33,1))," Dc=,",num2str(beta(35,1)),"±", num2str(betaBoot(35,1))," FD=,",num2str(beta(37,1)),"±", num2str(betaBoot(37,1)), "; Prob:",num2str(locModelProb)]);
        case 12
            yval = NAA(beta([38,40,39,41,42],1), dr,ri, tau, dimVal, Ntracks);
            col = colors(12,:);
            paramErrors = betaBoot([38,40,39,41,42],1);
            locModelProb = modelProb(modelProb(:,1) == 12,2);
            exportParamNameArray{i,1} = 12;
            exportParamNameArray{i,2} = ["Anomal+Anomal D1";"Anomal+Anomal alpha1";"Anomal+Anomal D2";"Anomal+Anomal alpha2";"Anomal+Anomal FD"];
            exportParamArray{i,1} = beta([38,40,39,41,42],1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Anomal+Anomal Diffusion, D1=',num2str(beta(38,1)), "±", num2str(betaBoot(38,1))," A1=,",num2str(beta(40,1)),"±", num2str(betaBoot(40,1))," D2=,",num2str(beta(39,1)),"±", num2str(betaBoot(39,1))," A2=,",num2str(beta(41,1)),"±", num2str(betaBoot(41,1))," FD=,",num2str(beta(42,1)),"±", num2str(betaBoot(42,1)), "; Prob:",num2str(locModelProb)]);
        case 13
            yval = NAConf(beta([43,44,45,47],1), dr,ri, tau, dimVal, Ntracks);
            col = colors(13,:);
            paramErrors = betaBoot([43,44,45,47],1);
            locModelProb = modelProb(modelProb(:,1) == 13,2);
            exportParamNameArray{i,1} = 13;
            exportParamNameArray{i,2} = ["Anomal+Confined Da";"Anomal+Confined alpha";"Anomal+Confined Dc";"Anomal+Confined FD"];
            exportParamArray{i,1} = beta([43,44,45,47],1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Anomal+Confined Diffusion, Da=',num2str(beta(43,1)), "±", num2str(betaBoot(43,1))," A=,",num2str(beta(44,1)),"±", num2str(betaBoot(44,1))," Dc=,",num2str(beta(45,1)),"±", num2str(betaBoot(45,1))," FD=,",num2str(beta(47,1)),"±", num2str(betaBoot(47,1)), "; Prob:",num2str(locModelProb)]);
        case 14
            yval = NConfConf(beta([48,50,52],1), dr,ri, tau, dimVal, Ntracks);
            col = colors(14,:);
            paramErrors = betaBoot([48,50,52],1);
            locModelProb = modelProb(modelProb(:,1) == 14,2);
            exportParamNameArray{i,1} = 14;
            exportParamNameArray{i,2} = ["Confined+Confined D1";"Confined+Confined D2";"Confined+Confined FD"];
            exportParamArray{i,1} = beta([48,50,52],1);
            exportParamArray{i,2} = paramErrors;
            exportParamArray{i,3} = locModelProb;
            legendLabelArray(i,1) = join(['Fit Confined+Confined Diffusion, D1=',num2str(beta(48,1)), "±", num2str(betaBoot(48,1))," D2=,",num2str(beta(50,1)),"±", num2str(betaBoot(50,1))," FD=,",num2str(beta(52,1)),"±", num2str(betaBoot(52,1)), "; Prob:",num2str(locModelProb)]);
    end
    %% do the plotting
    plot(Axes,ri,yval,"Color",col,'LineWidth',linWidth);

end
%% set the legend
legendLabelArray = ["Jump Distances";legendLabelArray];
legend(Axes,legendLabelArray.');

end



%% Extra functions for plotting
function [out] = ND(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
D = params;
dr = dr*Ntracks;
switch dim
    case 1
        out = dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau));
    case 2
        out = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
    case 3
        out = dr*ri.^2/(2*sqrt(pi)*(D*tau)^(3/2)).*exp(-ri.^2/(4*D*tau));
end
end

function [out] = NConf(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
dr = dr*Ntracks;
D = params(1);
maRad = params(2);
switch dim
    case 1
        %k = dr, -(erf((sqrt(d*t)*(2*r+k))/4))/(d*t)+(erf((sqrt(d*t)*(2*r-k))/4))/(d*t)+2*k
        %out = 0.5*(erf((-ri+dr/2)/(sqrt(4*D*tau)))-erf((-ri-dr/2)/(sqrt(4*D*tau))));
        %out = dr-((dr/sqrt(pi*D*tau))*exp((-ri.^2)/4*D*tau));
        % S = 0;
        % tol = 1e-4;
        % term = inf;
        % n = 1;
        % while abs(term) > tol
        %     term = exp(-((n*pi/2)^2)+tau*D/maRad^2) * cos(n*pi*ri/maRad*2) * cos(n*pi*0/maRad*2);
        %     S = S+term;
        %     n = n+1;
        % end
        % predicted2 = S;
        %out = dr/maRad*2+2*dr/maRad*2*predicted2;
        out = (dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau)));
    case 2
        %out = 0.5*(exp(-(((ri-dr/2).^2)/( 4*D*tau)))-exp(-(((ri+dr/2).^2)/( 4*D*tau))));
        %  k = dr, 2*pi*k*r-(erf(pi*k*sqrt(d*t)*r))/(4*sqrt(pi)*d*t*sqrt(d*t))
        out = (ri.*dr*2*pi) - ((ri.*dr)/(2*D*tau)).*exp((-ri.^2)/(4*D*tau));
    case 3
        %out =(ri/sqrt(4*D*tau)).*(exp(-(((ri-dr/2).^2)/( 4*D*tau)))-exp(-(((ri+dr/2).^2)/( 4*D*tau))));
        %-(sqrtd*sqrt(1/(d*t))*sqrt(t)*erf((sqrt(1/(d*t))*r)/2))/(4*pi)+r/(4*pi*sqrt(pi)*sqrtd*sqrt(t)*e^((r^2)/(4*d*t)))+(r^3)/3
        out = (4*pi*ri.^2*dr) - ((ri.^2)*dr/(2*sqrt(pi)*(D*tau)^(3/2))).*exp(-ri.^2/(4*D*tau));
end
end

function [out] = NV(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
dr = dr*Ntracks;
%Input Parameters
V = params(1);
D = params(2);

switch dim
    case 1
        z = -(ri.^2+V^2*tau^2)/(4*D*tau);
        y = ri*V/(2*D);
        out =dr/((4*pi*D*tau)^(1/2)).*exp(z+y)+dr/((4*pi*D*tau)^(1/2)).*exp(z-y);
    case 2
        z = -(ri.^2+V^2*tau^2)/(4*D*tau);
        y = ri*V/(2*D);
        out = dr*ri/(2*D*tau).*exp(z).*besseli(0, y);
    case 3
        z = -(ri.^2+V^2*tau^2)/(4*D*tau);
        y = ri*V/(2*D);
        out =dr*ri.^2*4*pi/((4*pi*D*tau)^(3/2)).*exp(z).*sinh(y)./y;
end
end

function [out] = NA(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
dr = dr*Ntracks;
%Input Parameters
Dalpha = params(1);
alpha = params(2);

%Setting integration limits based on alpha


switch dim
    case 1
        if alpha < 0.5
            min=-300^(.5/alpha); %limits on inverse laplace transform
        else
            min=-500;
        end

        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*...
            exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
        out = dr/((Dalpha)^(1/2)).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-5,'RelTol',1e-3));

    case 2
        if alpha < 0.5
            min=-300^(.5/alpha); %limits on inverse laplace transform
        else
            min=-500;
        end
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        out=dr*ri/(Dalpha).*abs(integral(fun,min,-1*min,...
            'ArrayValued',true,'AbsTol',1e-6));
    case 3
        if alpha < 0.5 && alpha>0.4
            min=-700; %limits on inverse laplace transform
        elseif alpha<= 0.4 && alpha>0.3
            min=-900;
        elseif alpha<=0.3
            min=-1100;
        else
            min=-700;
        end
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        intsize=200;
        numint=ceil(abs(min+100)/intsize);
        predicted=zeros(size(ri));
        for i=1:numint
            predicted=predicted+2*dr*ri/(Dalpha).*abs(integral(fun,min+(i-1)...
                *intsize-1i*1e-6,min+i*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        out=predicted+dr*ri/(Dalpha).*abs(integral(fun,-100-1i*1e-6,...
            100-1i*1e-6,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
end
end

function [out] = NDD(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
dr = dr*Ntracks;
%Input Parameter
D1=params(1);
D2=params(2);
fd=params(3);

switch dim
    case 1
        predicted1=dr/((pi*D1*tau)^(1/2)).*exp(-ri.^2/(4*D1*tau));
        predicted2=dr/((pi*D2*tau)^(1/2)).*exp(-ri.^2/(4*D2*tau));
        out=fd*predicted1+(1-fd)*predicted2;
    case 2
        predicted1 = dr*ri/(2*D1*tau).*exp(-ri.^2/(4*D1*tau));
        predicted2 = dr*ri/(2*D2*tau).*exp(-ri.^2/(4*D2*tau));
        out=fd*predicted1+(1-fd)*predicted2;
    case 3
        predicted1 = dr*ri.^2/(2*sqrt(pi)*(D1*tau)^(3/2)).*exp(-ri.^2/(4*D1*tau));
        predicted2 = dr*ri.^2/(2*sqrt(pi)*(D2*tau)^(3/2)).*exp(-ri.^2/(4*D2*tau));
        out=fd*predicted1+(1-fd)*predicted2;
end
end

function [out] = NDV(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
dr = dr*Ntracks;
%Input Parameters
D=params(1);
V=params(2);
Dv=params(3);
fd=params(4);

switch dim
    case 1
        z = -(ri.^2+V^2*tau^2)/(4*Dv*tau);
        y = ri*V/(2*Dv);
        out=fd*dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau))+...
            (1-fd)*(dr/((4*pi*Dv*tau)^(1/2)).*exp(z+y)+dr/((4*pi*Dv*tau)^(1/2)).*exp(z-y));
    case 2
        predicted1 = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
        z = -(ri.^2+V^2*tau^2)/(4*Dv*tau);
        y = ri*V/(2*Dv);
        predicted2 = dr*ri/(2*Dv*tau).*exp(z).*besseli(0, y);

        out=fd*predicted1+(1-fd)*predicted2;
    case 3
        predicted1 = dr*ri.^2/(2*sqrt(pi)*(D*tau)^(3/2)).*exp(-ri.^2/(4*D*tau));
        z = -(ri.^2+V^2*tau^2)/(4*Dv*tau);
        y = ri*V/(2*Dv);
        predicted2=dr*ri.^2*4*pi/((4*pi*Dv*tau)^(3/2)).*exp(z).*sinh(y)./y;


        out=fd*predicted1+(1-fd)*predicted2;
end
end

function [out] = NVV(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
dr = dr*Ntracks;
%Input Parameters
D1=params(1);
V1=params(2);
D2=params(3);
V2=params(4);
fd=params(5);

switch dim
    case 1
        z1 = -(ri.^2+V1^2*tau^2)/(4*D1*tau);
        y1 = ri*V1/(2*D1);
        z2 = -(ri.^2+V2^2*tau^2)/(4*D2*tau);
        y2 = ri*V2/(2*D2);
        out=fd*(dr/((4*pi*D1*tau)^(1/2)).*exp(z1+y1)+dr/((4*pi*D1*tau)^(1/2)).*exp(z1-y1))+...
            (1-fd)*(dr/((4*pi*D2*tau)^(1/2)).*exp(z2+y2)+dr/((4*pi*D2*tau)^(1/2)).*exp(z2-y2));
    case 2
        z1 = -(ri.^2+V1^2*tau^2)/(4*D1*tau);
        y1 = ri*V1/(2*D1);
        z2 = -(ri.^2+V2^2*tau^2)/(4*D2*tau);
        y2 = ri*V2/(2*D2);
        predicted1 = dr*ri/(2*D1*tau).*exp(z1).*besseli(0, y1);
        predicted2 = dr*ri/(2*D2*tau).*exp(z2).*besseli(0, y2);
        out=fd*predicted1+(1-fd)*predicted2;
    case 3
        z1 = -(ri.^2+V1^2*tau^2)/(4*D1*tau);
        y1 = ri*V1/(2*D1);
        predicted1=dr*ri.^2*4*pi/((4*pi*D1*tau)^(3/2)).*exp(z1).*sinh(y1)./y1;
        z2 = -(ri.^2+V2^2*tau^2)/(4*D2*tau);
        y2 = ri*V2/(2*D2);
        predicted2=dr*ri.^2*4*pi/((4*pi*D2*tau)^(3/2)).*exp(z2).*sinh(y2)./y2;
        out=fd*predicted1+(1-fd)*predicted2;
end

end

function [out] = NDA(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
dr = dr*Ntracks;
%Input Parameters
D=params(1);
Dalpha = params(2);
alpha = params(3);
fd=params(4);


switch dim
    case 1
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
        if alpha < 0.5
            min=-300^(.5/alpha);
        else
            min=-500;
        end
        %Plot the predicted JDD on top of it
        out=fd*dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau))+...
            (1-fd)*(dr/((Dalpha)^(1/2)).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',2e-5)));
    case 2
        if alpha < 0.5
            min1=-300^(.5/alpha);
        else
            min1=-500;
        end
        %Plot the predicted JDD on top of it
        predicted1 = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        predicted2=dr*ri/(Dalpha).*abs(integral(fun,min1,-1*min1,...
            'ArrayValued',true,'AbsTol',1e-6));

        out=fd*predicted1+(1-fd)*predicted2;
    case 3
        if alpha < 0.5 && alpha>0.4
            min=-700; %limits on inverse laplace transform
        elseif alpha<= 0.4 && alpha>0.3
            min=-900;
        elseif alpha<=0.3
            min=-1100;
        else
            min=-700;
        end

        %predicted JDD probabilities based on input parameters
        %For 3D anomalous diffusion, we found that splitting up the
        %integration into chunks led to the most accurate result.
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        intsize=200;
        numint=ceil(abs(min+100)/intsize);
        predicted=zeros(size(ri));
        for i=1:numint
            predicted=predicted+2*dr*ri/(Dalpha).*abs(integral(fun,min+(i-1)...
                *intsize-1i*1e-6,min+i*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        predicted2=predicted+dr*ri/(Dalpha).*abs(integral(fun,-100-1i*1e-6,...
            100-1i*1e-6,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));

        %Plot the predicted JDD on top of it
        predicted1 = dr*ri.^2/(2*sqrt(pi)*(D*tau)^(3/2)).*exp(-ri.^2/(4*D*tau));

        out=fd*predicted1+(1-fd)*predicted2;
end
end

function [out] = NDConf(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
dr = dr*Ntracks;
%Input Parameters
D=params(1);
DConf = params(2);
fd=params(3);




switch dim
    case 1
       % out=fd*dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau))+...
       %     (1-fd)*(0.5*(erf((-ri+dr/2)/(sqrt(4*DConf*tau)))-erf((-ri-dr/2)/(sqrt(4*DConf*tau)))));
       out=fd*dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau))+...
            (1-fd)*(dr-((dr/sqrt(pi*DConf*tau))*exp((-ri.^2)/4*DConf*tau)));
    case 2
        predicted1 = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
        %predicted2 = 0.5*(exp(-(((ri-dr/2).^2)/( 4*DConf*tau)))-exp(-(((ri+dr/2).^2)/( 4*DConf*tau))));
        predicted2 = (ri.*dr*2*pi) - ((ri.*dr)/(2*DConf*tau)).*exp((-ri.^2)/(4*DConf*tau));
        out=fd*predicted1+(1-fd)*predicted2;

    case 3
        predicted1 = dr*ri.^2/(2*sqrt(pi)*(D*tau)^(3/2)).*exp(-ri.^2/(4*D*tau));
       % predicted2 = (ri/sqrt(4*DConf*tau)).*(exp(-(((ri-dr/2).^2)/( 4*DConf*tau)))-exp(-(((ri+dr/2).^2)/( 4*DConf*tau))));
       predicted2 = (4*pi*ri.^2*dr) - ((ri.^2)*dr/(2*sqrt(pi)*(DConf*tau)^(3/2))).*exp(-ri.^2/(4*DConf*tau));

        out=fd*predicted1+(1-fd)*predicted2;
end
end

function [out] = NDirAno(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
%[msdDV,msdV,msdA,msdDA,.5]
%Input Parameters
DDir=params(1);
DV = params(2);
Dalpha = params(4);
alpha = params(3);
fd=params(5);
dr = dr*Ntracks;
switch dim
    case 1
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
        if alpha < 0.5
            min=-300^(.5/alpha);
        else
            min=-500;
        end
        %Plot the predicted JDD on top of it
        z1 = -(ri.^2+DV^2*tau^2)/(4*DDir*tau);
        y1 = ri*DV/(2*DDir);
        out=fd*(dr/((4*pi*DDir*tau)^(1/2)).*exp(z1+y1)+dr/((4*pi*DDir*tau)^(1/2)).*exp(z1-y1))+...
            (1-fd)*(dr/((Dalpha)^(1/2)).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',2e-5)));

    case 2
        if alpha < 0.5
            min=-300^(.5/alpha);
        else
            min=-500;
        end
        %Plot the predicted JDD on top of it
        z = -(ri.^2+DV^2*tau^2)/(4*DDir*tau);
        y = ri*DV/(2*DDir);
        predicted1 = dr*ri/(2*DDir*tau).*exp(z).*besseli(0, y);

        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        predicted2=dr*ri/(Dalpha).*abs(integral(fun,min,-1*min,...
            'ArrayValued',true,'AbsTol',1e-6));

        out=fd*predicted1+(1-fd)*predicted2;

    case 3
        z = -(ri.^2+DV^2*tau^2)/(4*DDir*tau);
        y = ri*DV/(2*DDir);
        predicted1=dr*ri.^2*4*pi/((4*pi*DDir*tau)^(3/2)).*exp(z).*sinh(y)./y;

        if alpha < 0.5 && alpha>0.4
            min=-700; %limits on inverse laplace transform
        elseif alpha<= 0.4 && alpha>0.3
            min=-900;
        elseif alpha<=0.3
            min=-1100;
        else
            min=-700;
        end

        %predicted JDD probabilities based on input parameters
        %For 3D anomalous diffusion, we found that splitting up the
        %integration into chunks led to the most accurate result.
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        intsize=200;
        numint=ceil(abs(min+100)/intsize);
        predicted=zeros(size(ri));
        for i=1:numint
            predicted=predicted+2*dr*ri/(Dalpha).*abs(integral(fun,min+(i-1)...
                *intsize-1i*1e-6,min+i*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        predicted2=predicted+dr*ri/(Dalpha).*abs(integral(fun,-100-1i*1e-6,...
            100-1i*1e-6,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));

        out=fd*predicted1+(1-fd)*predicted2;
end
end

function [out] = NDirConf(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
%[msdDV,msdV,msdA,msdDA,.5]
%Input Parameters
% x0 = [msdDV,msdV,msdConD,0.5];
DDir=params(1);
DV = params(2);
DConf = params(3);
fd=params(4);
dr = dr*Ntracks;







switch dim
    case 1
        z1 = -(ri.^2+DV^2*tau^2)/(4*DDir*tau);
        y1 = ri*DV/(2*DDir);

        %Plot the predicted JDD on top of it
        out=fd*(dr/((4*pi*DDir*tau)^(1/2)).*exp(z1+y1)+dr/((4*pi*DDir*tau)^(1/2)).*exp(z1-y1))+...
            (1-fd)*(dr-((dr/sqrt(pi*DConf*tau))*exp((-ri.^2)/4*DConf*tau)));
    case 2
        z = -(ri.^2+DV^2*tau^2)/(4*DDir*tau);
        y = ri*DV/(2*DDir);
        predicted1 = dr*ri/(2*DDir*tau).*exp(z).*besseli(0, y);
        %predicted2 = 0.5*(exp(-(((ri-dr/2).^2)/( 4*Dconf*tau)))-exp(-(((ri+dr/2).^2)/( 4*Dconf*tau))));
        predicted2 = (ri.*dr*2*pi) - ((ri.*dr)/(2*DConf*tau)).*exp((-ri.^2)/(4*DConf*tau));

        %Plot the predicted JDD on top of it
        out=fd*predicted1+(1-fd)*predicted2;

    case 3
        z = -(ri.^2+DV^2*tau^2)/(4*DDir*tau);
        y = ri*DV/(2*DDir);
        predicted1=dr*ri.^2*4*pi/((4*pi*DDir*tau)^(3/2)).*exp(z).*sinh(y)./y;
        predicted2 = (4*pi*ri.^2*dr) - ((ri.^2)*dr/(2*sqrt(pi)*(DConf*tau)^(3/2))).*exp(-ri.^2/(4*DConf*tau));

        %Plot the predicted JDD on top of it
        out=fd*predicted1+(1-fd)*predicted2;
end
end

function [out] = NAA(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
dr = dr*Ntracks;
%Input Parameters
Dalpha1 = params(1);
alpha1 = params(2);
Dalpha2 = params(3);
alpha2 = params(4);
fd=params(5);

switch dim
    case 1
        if alpha1 < 0.5
            min1=-300^(.5/alpha1); %limits on inverse laplace transform
        else
            min1=-500;
        end
        if alpha2 < 0.5
            min2=-300^(.5/alpha2); %limits on inverse laplace transform
        else
            min2=-500;
        end

        %predicted JDD probabilities based on input parameters
        fun1=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha1/2-1)/(2.*pi).*...
            exp(-ri./(sqrt(Dalpha1)).*((1i*p)^(alpha1/2)));
        fun2=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha2/2-1)/(2.*pi).*...
            exp(-ri./(sqrt(Dalpha2)).*((1i*p)^(alpha2/2)));

        out=fd*(dr/((Dalpha1)^(1/2)).*abs(integral(fun1,min1,-1*min1,'ArrayValued',true,'AbsTol',2e-5)))+...
            (1-fd)*(dr/((Dalpha2)^(1/2)).*abs(integral(fun2,min2,-1*min2,'ArrayValued',true,'AbsTol',2e-5)));
    case 2
        %Setting integration limits based on alpha
        if alpha1 < 0.5
            min1=-300^(.5/alpha1); %limits on inverse laplace transform
        else
            min1=-500;
        end
        if alpha2 < 0.5
            min2=-300^(.5/alpha2); %limits on inverse laplace transform
        else
            min2=-500;
        end

        %predicted JDD probabilities based on input parameters
        fun1=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha1-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(Dalpha1)).*((1i*p)^(alpha1/2))));
        predicted1=dr*ri/(Dalpha1).*abs(integral(fun1,min1,-1*min1,...
            'ArrayValued',true,'AbsTol',1e-6));

        fun2=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha2-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(Dalpha2)).*((1i*p)^(alpha2/2))));
        predicted2=dr*ri/(Dalpha2).*abs(integral(fun2,min2,-1*min2,...
            'ArrayValued',true,'AbsTol',1e-6));

        out=fd*predicted1+(1-fd)*predicted2;

    case 3
        if alpha1 < 0.5 && alpha1>0.4
            min1=-700; %limits on inverse laplace transform
        elseif alpha1<= 0.4 && alpha1>0.3
            min1=-900;
        elseif alpha1<=0.3
            min1=-1100;
        else
            min1=-700;
        end

        %predicted JDD probabilities based on input parameters
        %For 3D anomalous diffusion, we found that splitting up the
        %integration into chunks led to the most accurate result.
        fun1=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha1-1)/(2.*pi).*...
            (exp(-ri./(sqrt(Dalpha1)).*((1i*p)^(alpha1/2))));
        intsize=200;
        numint=ceil(abs(min1+100)/intsize);
        predicted=zeros(size(ri));
        for i=1:numint
            predicted=predicted+2*dr*ri/(Dalpha1).*abs(integral(fun1,min1+(i-1)...
                *intsize-1i*1e-6,min1+i*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        predicted1=predicted+dr*ri/(Dalpha1).*abs(integral(fun1,-100-1i*1e-6,...
            100-1i*1e-6,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));


        if alpha2 < 0.5 && alpha2>0.4
            min2=-700; %limits on inverse laplace transform
        elseif alpha2<= 0.4 && alpha2>0.3
            min2=-900;
        elseif alpha2<=0.3
            min2=-1100;
        else
            min2=-700;
        end

        %predicted JDD probabilities based on input parameters
        %For 3D anomalous diffusion, we found that splitting up the
        %integration into chunks led to the most accurate result.
        fun2=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha2-1)/(2.*pi).*...
            (exp(-ri./(sqrt(Dalpha2)).*((1i*p)^(alpha2/2))));
        intsize=200;
        numint=ceil(abs(min2+100)/intsize);
        predicted=zeros(size(ri));
        for i=1:numint
            predicted=predicted+2*dr*ri/(Dalpha2).*abs(integral(fun2,min2+(i-1)...
                *intsize-1i*1e-6,min2+i*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        predicted2=predicted+dr*ri/(Dalpha2).*abs(integral(fun2,-100-1i*1e-6,...
            100-1i*1e-6,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));

        out=fd*predicted1+(1-fd)*predicted2;
end
end

function [out] = NAConf(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
%x0 = [msdA,msdDA,msdConD,0.5];
%Input Parameters
Dalpha1 = params(1);
alpha1 = params(2);
DConf = params(3);
fd=params(4);
dr = dr*Ntracks;



switch dim
    case 1
        if alpha1 < 0.5
            min1=-300^(.5/alpha1); %limits on inverse laplace transform
        else
            min1=-500;
        end


        %predicted JDD probabilities based on input parameters
        fun1=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha1/2-1)/(2.*pi).*...
            exp(-ri./(sqrt(Dalpha1)).*((1i*p)^(alpha1/2)));

        out=fd*(dr/((Dalpha1)^(1/2)).*abs(integral(fun1,min1,-1*min1,'ArrayValued',true,'AbsTol',2e-5)))+...
            (1-fd)*(dr-((dr/sqrt(pi*DConf*tau))*exp((-ri.^2)/4*DConf*tau)));
    case 2
        %Setting integration limits based on alpha
        if alpha1 < 0.5
            min1=-300^(.5/alpha1); %limits on inverse laplace transform
        else
            min1=-500;
        end


        %predicted JDD probabilities based on input parameters
        fun1=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha1-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(Dalpha1)).*((1i*p)^(alpha1/2))));
        predicted1=dr*ri/(Dalpha1).*abs(integral(fun1,min1,-1*min1,...
            'ArrayValued',true,'AbsTol',1e-6));
        predicted2 = (ri.*dr*2*pi) - ((ri.*dr)/(2*DConf*tau)).*exp((-ri.^2)/(4*DConf*tau));

        out=fd*predicted1+(1-fd)*predicted2;

    case 3
        if alpha1 < 0.5 && alpha1>0.4
            min1=-700; %limits on inverse laplace transform
        elseif alpha1<= 0.4 && alpha1>0.3
            min1=-900;
        elseif alpha1<=0.3
            min1=-1100;
        else
            min1=-700;
        end

        %predicted JDD probabilities based on input parameters
        %For 3D anomalous diffusion, we found that splitting up the
        %integration into chunks led to the most accurate result.
        fun1=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha1-1)/(2.*pi).*...
            (exp(-ri./(sqrt(Dalpha1)).*((1i*p)^(alpha1/2))));
        intsize=200;
        numint=ceil(abs(min1+100)/intsize);
        predicted=zeros(size(ri));
        for i=1:numint
            predicted=predicted+2*dr*ri/(Dalpha1).*abs(integral(fun1,min1+(i-1)...
                *intsize-1i*1e-6,min1+i*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        predicted1=predicted+dr*ri/(Dalpha1).*abs(integral(fun1,-100-1i*1e-6,...
            100-1i*1e-6,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));



        predicted2 = (4*pi*ri.^2*dr) - ((ri.^2)*dr/(2*sqrt(pi)*(DConf*tau)^(3/2))).*exp(-ri.^2/(4*DConf*tau));

        out=fd*predicted1+(1-fd)*predicted2;
end

end

function [out] = NConfConf(params, dr,ri, tau, dim, Ntracks)
%Generates the values needed to perform the LSQ non-linear fit for
%the directed motion model against a given JDD.
%[msdConD*0.9,msdConD*1.1,0.5]
%Input Parameter
D1=params(1);
D2=params(2);
fd=params(3);
dr = dr*Ntracks;
switch dim
    case 1
        out=fd*(dr-((dr/sqrt(pi*D1*tau))*exp((-ri.^2)/4*D1*tau)))+( ...
            1-fd)*(dr-((dr/sqrt(pi*D2*tau))*exp((-ri.^2)/4*D2*tau)));
    case 2
        predicted1 = (ri.*dr*2*pi) - ((ri.*dr)/(2*D1*tau)).*exp((-ri.^2)/(4*D1*tau));
        predicted2 = (ri.*dr*2*pi) - ((ri.*dr)/(2*D2*tau)).*exp((-ri.^2)/(4*D2*tau));
        out=fd*predicted1+(1-fd)*predicted2;

    case 3
        predicted1 = (4*pi*ri.^2*dr) - ((ri.^2)*dr/(2*sqrt(pi)*(D1*tau)^(3/2))).*exp(-ri.^2/(4*D1*tau));
        predicted2 = (4*pi*ri.^2*dr) - ((ri.^2)*dr/(2*sqrt(pi)*(D2*tau)^(3/2))).*exp(-ri.^2/(4*D2*tau));
        out=fd*predicted1+(1-fd)*predicted2;
end

end

