%Parameter fitting for 1D Jump Distance Distributions
%Thomas Müller
%Last Updated: 1/8/24

%This function takes JDD data as the input and outputs the parameters for
%three different models using lsqcurvefit with weighting.
%our model fitting method is only valid for alpha <1. In the case that
%alpha >1, you can ignore fits, or deal with fitting in the model fitting
%step. As an alternative, you can impose bound constraints on any fits with
%anomalous diffusion.
%This function also allows for combination models, such as double
%diffusion, increasing computation time, but giving more options for fits.

%%%%%%%%%%INPUTS%%%%%%%%%%
%tau--the duration of each trajectory (the time lag*dt)
%dr--the width of each bin in the JDD histogram
%ri--vector of the midpoints of the JDD bins
%yi--vector of the percentage of trajectories in each JDD bin
%N--the number of trajectories
%points--the number of time points in each trajectory
%dt--the time step
%x1--the 1D trajectory

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%param: strut containing parameters determined by the LSQ fit for each model.

function param = ExtendedModelFitting2DwithComboModels(modelsToTest, windowType, tau, dr, ri, yi, Ni ,N, points, dt, x1, x2)
%Set up optimization options: can be customized for better fits (espeically
%tolerances)
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'MaxFunctionEvaluations',1000,'FunctionTolerance',1e-4,...
    'StepTolerance',1e-4,'Display','off');

%Create struct to hold all the parameters
param = struct('FreeD',NaN,'DirV', NaN, 'DirDv', NaN,'AnoDalpha',NaN,'Anoalpha',NaN,...
    'FFD1', NaN,'FFD2',NaN,'FFfdDD',NaN,...
    'ConfD', NaN', 'ConfRad', NaN', 'FreeDirFreeD', NaN,'FreeDirV', NaN,'FreeDirVD', NaN,'FreeDirfdDV', NaN,...
    'FreeAnomFreeD', NaN,'FreeAnomAnomD', NaN,'FreeAnomalphada', NaN,'FreeAnomfdDA', NaN, ...
    'FreeConfFD', NaN, 'FreeConfConD', NaN,'FreeConfConR', NaN,'FreeConfConfd', NaN, ...
    'DirDirD1', NaN,'DirDirD2', NaN,'DirDirV1', NaN,'DirDirV2', NaN,'DirDirfd', NaN, ...
    'DirAnomD1', NaN,'DirAnomV1', NaN,'DirAnomD2', NaN,'DirAnomAlph2', NaN,'DirAnomfd', NaN, ...
    'DirConfD1', NaN,'DirConfV1', NaN,'DirConfD2', NaN,'DirConfR2', NaN,'DirConffd', NaN, ...
    'AnomAnomD1', NaN,'AnomAnomD2', NaN,'AnomAnomA1', NaN,'AnomAnomA2', NaN,'AnomAnomfd', NaN, ...
    'AnomConfD1', NaN,'AnomConfA1', NaN,'AnomConfD2', NaN,'AnomConfR2', NaN,'AnomConffd', NaN,'ConfConfD1', ...
    NaN,'ConfConfR1', NaN,'ConfConfD2',NaN,'ConfConfR2', NaN,'ConfConffd', NaN);


%Weighting vector
%We chose to weight by 1/bin count probabilties (Ni/N=bin count
%probabilities). Each bin was given an artificial additional count so no
%bin was empty
wi=1./((Ni+1)/(N+length(Ni)));
weights=wi;
%if you are using a sliding construction for the jdd then you do NOT want
%weights.
if windowType
    weights=1; %for a correlated or sliding JDD.
end

%Getting initial seeds by using MSD Analysis
%make time vector
t=dt*(0:points);

%finding avg(x^2). Have to account for possible missing data. Outliers can
%mess this calculation up, but there isn't too much that can be done in
%that case, since using nanmedian causes other issues.
xsqavg=nanmean((x1-repmat(x1(1,:),size(x1,1),1)).^2+...
    (x2-repmat(x2(1,:),size(x2,1),1)).^2,2);

t=t';

for i = 1:size(modelsToTest,1)
    switch modelsToTest(i)
        case 1
            %1 = FreeDiff
            %Pure Diffusion Seeds
            [p]=polyfitZero(t, xsqavg,1);
            msdD=abs(p(1)/4);

            %Testing 1D Diffusion Model
            x0=msdD;
            [temp]=lsqnonlin(@ND,x0,[],[],options);
            param.FreeD=temp(1);
        case 2
            %2 = Directed
            %Directed Diffusion Seeding
            options2=optimset('Display','off');
            [p]=polyfitZero(t,xsqavg,2);
            if p(1)< 0
                paramF = fmincon(@(x) norm(x(1)^2*t.^2 + 4*x(2)*t - xsqavg),[0.1 0.1],...
                    [],[],[],[],[0 0],[],[],options2);
                msdV = paramF(1);
                msdDV = paramF(2);
            else
                msdV = sqrt(abs(p(1)));
                msdDV = abs(p(2)/4);
            end

            %Testing 1D Directed Motion Model
            x0 = [msdV, msdDV];
            [temp] = lsqnonlin(@NV, x0, [], [], options);
            param.DirV = temp(1);
            param.DirDv = temp(2);
        case 3
            %3 = Anomalous

            %Anomalous Diffusion Seeding
            [p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
            msdA=abs(p(1));
            msdDA=abs(exp(p(2)-log(4)+log(gamma(1+msdA))));

            %Testing 1D Anomalous Diffusion Model
            %we deal with any anomalous problems (i.e. alpha>1 in model fitting)
            x0 = [msdDA, msdA];
            %if you want to fully constrain alpha, you do it as below
            %[temp]=lsqnonlin(@NA,x0,[0 0],[500 1],options);
            %this will change the method used for curve fitting to the trust-region
            %reflective method, and can add computation time. We prefer to just exclude
            %models that fit over alpha=1.
            [temp]=lsqnonlin(@NA,x0,[],[],options);
            param.AnoDalpha=temp(1);
            param.Anoalpha=temp(2);
        case 4
            %4 = Confined
            %Confied diffusion seeding
            maxLengthx = max(abs(x1(end,:)-x1(1,:)));
            maxLengthy = max(abs(x2(end,:)-x2(1,:)));
            maxLength = max([maxLengthy, maxLengthx]);
            maxRad = maxLength/2;
            [p]=polyfitZero(t, xsqavg,1);
            msdD=abs(p(1)/2);
            x0 = [msdD, maxRad];
            %model fitting
            [temp]=lsqnonlin(@NConf,x0,[1,1],[Inf,maxLength],options);
            param.ConfD=temp(1);
            param.ConfRad=temp(2);
        case 5
            %5 = Free+Free
            %Pure Diffusion Seeds
             %Pure Diffusion Seeds
            [p]=polyfitZero(t, xsqavg,1);
            msdD=abs(p(1)/4);

            %Testing 1D Double Diffusion Model
            x0=[0.5*msdD,1.5*msdD,.5]; %so they start out different
            [temp]=lsqnonlin(@NDD,x0,[],[],options);
            param.FFD1=temp(1);
            param.FFD2=temp(2);
            param.FFfdDD=temp(3);
        case 6
            %6 = Free+Directed
             %Pure Diffusion Seeds
            [p]=polyfitZero(t, xsqavg,1);
            msdD=abs(p(1)/4);
             %Directed Diffusion Seeding
            options2=optimset('Display','off');
            [p]=polyfitZero(t,xsqavg,2);
            if p(1)< 0
                paramF = fmincon(@(x) norm(x(1)^2*t.^2 + 4*x(2)*t - xsqavg),[0.1 0.1],...
                    [],[],[],[],[0 0],[],[],options2);
                msdV = paramF(1);
                msdDV = paramF(2);
            else
                msdV = sqrt(abs(p(1)));
                msdDV = abs(p(2)/4);
            end
            %Testing 1D Diffusion Directed Combo Model
            x0=[msdD,msdV, msdDV,.5]; %so they start out different
            [temp]=lsqnonlin(@NDV,x0,[],[],options);
            param.FreeDirFreeD=temp(1);
            param.FreeDirV=temp(2);
            param.FreeDirVD=temp(3);
            param.FreeDirfdDV=temp(4);
        case 7
            %7 = Free+Anomal
             %Pure Diffusion Seeds
            [p]=polyfitZero(t, xsqavg,1);
            msdD=abs(p(1)/4);

             %Anomalous Diffusion Seeding
            [p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
            msdA=abs(p(1));
            msdDA=abs(exp(p(2)-log(4)+log(gamma(1+msdA))));

            %Testing 1D Diffusion-Anomalous Combo Model
            x0=[msdD,msdDA,msdA,.5]; %so they start out different
            [temp]=lsqnonlin(@NDA,x0,[],[],options);
            param.FreeAnomFreeD=temp(1);
            param.FreeAnomAnomD=temp(2);
            param.FreeAnomalphada=temp(3);
            param.FreeAnomfdDA=temp(4);
        case 8
            %8 = Free+Confined
             %Pure Diffusion Seeds
            [p]=polyfitZero(t, xsqavg,1);
            msdD=abs(p(1)/4);
            %Confied diffusion seeding
            [p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
            msdConR=0.5; %confiend alpha
            msdConD=abs(exp(p(2)-log(4)+log(gamma(1+msdConR))));
            x0 = [msdD,msdConD,0.5];

            [temp]=lsqnonlin(@NDConf,x0,[],[],options);

            param.FreeConfFD=temp(1);
            param.FreeConfConD=temp(2);
            param.FreeConfConR=temp(2)*tau;
            param.FreeConfConfd=temp(3);

        case 9
            %9 = Directed+Directed
            %Directed Diffusion Seeding
            options2=optimset('Display','off');
            [p]=polyfitZero(t,xsqavg,2);
            if p(1)< 0
                paramF = fmincon(@(x) norm(x(1)^2*t.^2 + 4*x(2)*t - xsqavg),[0.1 0.1],...
                    [],[],[],[],[0 0],[],[],options2);
                msdV = paramF(1);
                msdDV = paramF(2);
            else
                msdV = sqrt(abs(p(1)));
                msdDV = abs(p(2)/4);
            end

            %Testing 1D Directed Motion Model
            x0 = [msdDV*0.9, msdV*0.9, msdDV*1.1, msdV*1.1, 0.5];
            [temp] = lsqnonlin(@NVV, x0, [], [], options);
            param.DirDirD1 = temp(1);
            param.DirDirV1 = temp(2);
            param.DirDirD2 = temp(3);
            param.DirDirV2 = temp(4);
            param.DirDirfd = temp(5);

        case 10
            %10 = Directed+Anomal
             %Directed Diffusion Seeding
            options2=optimset('Display','off');
            [p]=polyfitZero(t,xsqavg,2);
            if p(1)< 0
                paramF = fmincon(@(x) norm(x(1)^2*t.^2 + 4*x(2)*t - xsqavg),[0.1 0.1],...
                    [],[],[],[],[0 0],[],[],options2);
                msdV = paramF(1);
                msdDV = paramF(2);
            else
                msdV = sqrt(abs(p(1)));
                msdDV = abs(p(2)/4);
            end
             %Anomalous Diffusion Seeding
            [p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
            msdA=abs(p(1));
            msdDA=abs(exp(p(2)-log(4)+log(gamma(1+msdA))));


            x0=[msdDV,msdV,msdA,msdDA,.5]; %so they start out different
            [temp]=lsqnonlin(@NDirAno,x0,[],[],options);
            param.DirAnomD1=temp(1);
            param.DirAnomV1=temp(2);
            param.DirAnomAlph2=temp(3);
            param.DirAnomD2=temp(4);
            param.DirAnomfd=temp(5);
        case 11
            %11 = Directed+Confined
            %Directed Diffusion Seeding
            options2=optimset('Display','off');
            [p]=polyfitZero(t,xsqavg,2);
            if p(1)< 0
                paramF = fmincon(@(x) norm(x(1)^2*t.^2 + 4*x(2)*t - xsqavg),[0.1 0.1],...
                    [],[],[],[],[0 0],[],[],options2);
                msdV = paramF(1);
                msdDV = paramF(2);
            else
                msdV = sqrt(abs(p(1)));
                msdDV = abs(p(2)/4);
            end
            %Confied diffusion seeding
            [p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
            msdConR=0.5; %confiend alpha
            msdConD=abs(exp(p(2)-log(4)+log(gamma(1+msdConR))));
            x0 = [msdDV,msdV,msdConD,0.5];
            [temp]=lsqnonlin(@NDirConf,x0,[],[],options);
            param.DirConfD1=temp(1);
            param.DirConfV1=temp(2);
            param.DirConfD2=temp(3);
            param.DirConfR2=temp(3)*tau;
            param.DirConffd=temp(4);
        case 12
            %12 = Anomal+Anomal
             %Anomalous Diffusion Seeding
            [p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
            msdA=abs(p(1));
            msdDA=abs(exp(p(2)-log(4)+log(gamma(1+msdA))));

            %Testing 1D Anomalous Diffusion Model
            %we deal with any anomalous problems (i.e. alpha>1 in model fitting)
            x0 = [msdDA*0.9, msdA*0.9,msdDA*1.1, msdA*1.1,0.5];
            %if you want to fully constrain alpha, you do it as below
            %[temp]=lsqnonlin(@NA,x0,[0 0],[500 1],options);
            %this will change the method used for curve fitting to the trust-region
            %reflective method, and can add computation time. We prefer to just exclude
            %models that fit over alpha=1.
            [temp]=lsqnonlin(@NAA,x0,[],[],options);
            param.AnomAnomD1=temp(1);
            param.AnomAnomA1=temp(2);
            param.AnomAnomD2=temp(3);
            param.AnomAnomA2=temp(4);
            param.AnomAnomfd=temp(5);
        case 13
            %13 = Anomal+Confined
             %Anomalous Diffusion Seeding
            [p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
            msdA=abs(p(1));
            msdDA=abs(exp(p(2)-log(4)+log(gamma(1+msdA))));

            %Confied diffusion seeding
            [p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
            msdConR=0.5; %confiend alpha
            msdConD=abs(exp(p(2)-log(4)+log(gamma(1+msdConR))));
            x0 = [msdDA,msdA,msdConD,0.5];
            [temp]=lsqnonlin(@NAConf,x0,[],[],options);
            param.AnomConfA1=temp(1);
            param.AnomConfD1=temp(2);
            param.AnomConfD2=temp(3);
            param.AnomConfR2=temp(3)*tau;
            param.AnomConffd=temp(4);
        case 14
            %14 = Confined+Confined
            %Confied diffusion seeding
            [p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
            msdConR=0.5; %confiend alpha
            msdConD=abs(exp(p(2)-log(4)+log(gamma(1+msdConR))));
            x0 = [msdConD*0.9,msdConD*1.1,0.5];
            [temp]=lsqnonlin(@NConfConf,x0,[],[],options);
            param.ConfConfD1=temp(1);
            param.ConfConfR1=temp(1)*tau;
            param.ConfConfD2=temp(2);
            param.ConfConfR2=temp(2)*tau;
            param.ConfConffd=temp(3);
    end
end


%% Fitting functions for the different modes
    function [outvals] = ND(x) %done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameter
        D=x(1);

        %predicted JDD probabilities based on input parameters
        predicted = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));

        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NConf(x) %done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameter
        D=x(1);
        maRad = x(2);

        %predicted JDD probabilities based on input parameters
        %predicted=(ri.*dr*2*pi)-((ri.*dr)/(D*tau)).*exp(-ri.^2/(4*D*tau));
        %predicted = (ri.*dr*2*pi) - ((ri.*dr)/(2*D*tau)).*exp((-ri.^2)/(4*D*tau));
        %predicted1 = dr/maRad^2+dr/maRad^2;
        
        predicted = conditional_probability_no_angles_confined(ri, 0.000001, tau, maRad, D);

        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end
    
    function P = conditional_probability_no_angles_confined(r, r_prime, t, a, D, N, M)
        % Parameters
        if nargin < 7
            M = 100; % Number of terms in the radial summation
        end
        if nargin < 6
            N = 100; % Number of terms in the angular summation
        end

        % Time scale
        tau1 = a^2 / D;

        % Initialize probability
        P = dr / (pi * a^2);

        % Summation over n
        for n = -N:N
            % Summation over m
            for m = 1:M
                alpha_nm = find_besselj_zero_derivative(n);

                radial_part = (alpha_nm^2 / (alpha_nm^2 - n^2)) * ...
                    exp(-alpha_nm^2 * t / tau1) * ...
                    (besselj(n, alpha_nm * r / a) * besselj(n, alpha_nm * r_prime / a) / ...
                    (besselj(n, alpha_nm)^2));

                P = P + (dr / (pi * a^2)) * radial_part;
            end
        end
    end

    function zero = find_besselj_zero_derivative(n)
        % Systematically find the n-th positive zero of the derivative of the Bessel function J_n
        % Initial guess
        x0 = (n - 0.5) * pi;
        delta = 0.1; % Smaller step size for searching interval
        max_iter = 10000; % Maximum iterations for searching

        % Initialize the search interval
        x1 = x0 - delta;
        x2 = x0 + delta;

        for i = 1:max_iter
            f1 = besselj(1, x1);
            f2 = besselj(1, x2);

            if f1 * f2 < 0
                % Interval [x1, x2] brackets a root
                zero = fzero(@(x) besselj(1, x), [x1, x2]);
                return;
            else
                % Move to the next interval
                x1 = x2;
                x2 = x2 + delta;
            end
        end

        % If no zero is found, return NaN
        zero = NaN;
    end



    function [outvals] = NV(x) %done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameters
        V = x(1);
        D = x(2);

        %predicted JDD probabilities based on input parameters
        z = -(ri.^2+V^2*tau^2)/(4*D*tau);
        y = ri*V/(2*D);
        predicted = dr*ri/(2*D*tau).*exp(z).*besseli(0, y);

        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NA(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameters
        Dalpha = x(1);
        alpha = x(2);

        %Setting integration limits based on alpha
        if alpha < 0.5
            min=-300^(.5/alpha); %limits on inverse laplace transform
        else
            min=-500;
        end

        %predicted JDD probabilities based on input parameters
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        predicted=dr*ri/(Dalpha).*abs(integral(fun,min,-1*min,...
            'ArrayValued',true,'AbsTol',1e-6));
        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NDD(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameter
        D1=x(1);
        D2=x(2);
        fd=x(3);

        %predicted JDD probabilities based on input parameters
        predicted1 = dr*ri/(2*D1*tau).*exp(-ri.^2/(4*D1*tau));
        predicted2 = dr*ri/(2*D2*tau).*exp(-ri.^2/(4*D2*tau));
        predicted=fd*predicted1+(1-fd)*predicted2;
        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NDV(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameters
        D=x(1);
        V=x(2);
        Dv=x(3);
        fd=x(4);

        %predicted JDD probabilities based on input parameters
        predicted1 = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
        z = -(ri.^2+V^2*tau^2)/(4*Dv*tau);
        y = ri*V/(2*Dv);
        predicted2 = dr*ri/(2*Dv*tau).*exp(z).*besseli(0, y);

        predicted=fd*predicted1+(1-fd)*predicted2;

        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NVV(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameters
        D1=x(1);
        V1=x(2);
        D2=x(3);
        V2=x(4);
        fd=x(5);

        %predicted JDD probabilities based on input parameters
        z1 = -(ri.^2+V1^2*tau^2)/(4*D1*tau);
        y1 = ri*V1/(2*D1);
        z2 = -(ri.^2+V2^2*tau^2)/(4*D2*tau);
        y2 = ri*V2/(2*D2);
        predicted1 = dr*ri/(2*D1*tau).*exp(z1).*besseli(0, y1);
        predicted2 = dr*ri/(2*D2*tau).*exp(z2).*besseli(0, y2);
        predicted=fd*predicted1+(1-fd)*predicted2;

        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NDA(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameters
        D=x(1);
        Dalpha = x(2);
        alpha = x(3);
        fd=x(4);

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

        predicted=fd*predicted1+(1-fd)*predicted2;

        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NDConf(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameters
        D=x(1);
        DConf = x(2);
        fd=x(3);

        %Plot the predicted JDD on top of it
        predicted1 = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
        predicted2 = (ri.*dr*2*pi) - ((ri.*dr)/(2*DConf*tau)).*exp((-ri.^2)/(4*DConf*tau));
        predicted=fd*predicted1+(1-fd)*predicted2;

        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NDirAno(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        %[msdDV,msdV,msdA,msdDA,.5]
        %Input Parameters
        DDir=x(1);
        DV = x(2);
        Dalpha = x(4);
        alpha = x(3);
        fd=x(5);

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

        predicted=fd*predicted1+(1-fd)*predicted2;


        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NDirConf(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        %[msdDV,msdV,msdA,msdDA,.5]
        %Input Parameters
        % x0 = [msdDV,msdV,msdConD,0.5];
        DDir=x(1);
        DV = x(2);
        DConf = x(3);
        fd=x(4);

        z = -(ri.^2+DV^2*tau^2)/(4*DDir*tau);
        y = ri*DV/(2*DDir);
        predicted1 = dr*ri/(2*DDir*tau).*exp(z).*besseli(0, y);
        predicted2 = (ri.*dr*2*pi) - ((ri.*dr)/(2*DConf*tau)).*exp((-ri.^2)/(4*DConf*tau));

        %Plot the predicted JDD on top of it
        predicted=fd*predicted1+(1-fd)*predicted2;
        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NAA(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.

        %Input Parameters
        Dalpha1 = x(1);
        alpha1 = x(2);
        Dalpha2 = x(3);
        alpha2 = x(4);
        fd=x(5);

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

        predicted=fd*predicted1+(1-fd)*predicted2;
        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NAConf(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        %x0 = [msdA,msdDA,msdConD,0.5];
        %Input Parameters
        Dalpha1 = x(1);
        alpha1 = x(2);
        DConf = x(3);
        fd=x(4);

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

        predicted=fd*predicted1+(1-fd)*predicted2;
        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NConfConf(x)%done
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        %[msdConD*0.9,msdConD*1.1,0.5]
        %Input Parameter
        D1=x(1);
        D2=x(2);
        fd=x(3);

        %predicted JDD probabilities based on input parameters
        predicted1 = (ri.*dr*2*pi) - ((ri.*dr)/(2*D1*tau)).*exp((-ri.^2)/(4*D1*tau));
        predicted2 = (ri.*dr*2*pi) - ((ri.*dr)/(2*D2*tau)).*exp((-ri.^2)/(4*D2*tau));
        predicted=fd*predicted1+(1-fd)*predicted2;
        %actual JDD probabilities
        actual = yi;

        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end
end


