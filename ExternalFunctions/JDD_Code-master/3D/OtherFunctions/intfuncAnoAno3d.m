%Integration function for 1D Pure Diffusion Directed Motion Combination
%Motion
%Rebecca Menssen
%Last Updated: 3/20/19

%This function creates the grid necessary to integrate across a range of
%diffusion parameters using trapz

%%%%%%%%%%INPUTS%%%%%%%%%%
%x1,x2,x3--equally spaced array of D1,D2, and fd values to integrate over
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%out--structure of probability values for a range of D values

function[out]=intfuncAnoAno3d(x1,x2,x3,x4, x5, N,yi,ri,dr,tau)
out=zeros(length(x1),length(x2),length(x3),length(x4),length(x5));
for i=1:length(x1)
    Dalpha1=x1(i);
    for k=1:length(x2)
        alpha1=x2(k);
        if alpha1 < 0.5 && alpha1>0.4
            min1=-700; %limits on inverse laplace transform
        elseif alpha1<= 0.4 && alpha1>0.3
            min1=-900;
        elseif alpha1<=0.3
            min1=-1100;
        else
            min1=-700;
        end

        %calculate predicted probabilities
        fun1=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha1-1)/(2.*pi).*...
            (exp(-ri./(sqrt(Dalpha1)).*((1i*p)^(alpha1/2))));
        intsize=200;
        numint=ceil(abs(min1+100)/intsize);
        zt1=zeros(size(ri));
        for o=1:numint
            zt1=zt1+2*dr*ri/(Dalpha1).*abs(integral(fun1,min1+(o-1)*intsize-1i*1e-6,min1+o*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        z1=zt1+dr*ri/(Dalpha1).*abs(integral(fun1,-100-1i*1e-6,100-1i*1e-6,...
            'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        for p=1:length(x3)
            Dalpha2=x3(p);
            for o=1:length(x4)
                alpha2=x4(o);
                if alpha2 < 0.5 && alpha2>0.4
                    min2=-700; %limits on inverse laplace transform
                elseif alpha2<= 0.4 && alpha2>0.3
                    min2=-900;
                elseif alpha2<=0.3
                    min2=-1100;
                else
                    min2=-700;
                end

                %calculate predicted probabilities
                fun2=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha2-1)/(2.*pi).*...
                    (exp(-ri./(sqrt(Dalpha2)).*((1i*p)^(alpha2/2))));
                intsize=200;
                numint=ceil(abs(min2+100)/intsize);
                zt2=zeros(size(ri));
                for q=1:numint
                    zt2=zt2+2*dr*ri/(Dalpha2).*abs(integral(fun2,min2+(q-1)*intsize-1i*1e-6,min2+q*(intsize)-1i*1e-6,...
                        'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
                end
                z2=zt2+dr*ri/(Dalpha2).*abs(integral(fun2,-100-1i*1e-6,100-1i*1e-6,...
                    'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));

                parfor l=1:length(x4)
                    [i,k,p,o,l]
                    fd=x4(l);
                    z=fd*z1+(1-fd)*z2;
                    %calculate P(JDD|model,parameters)
                    denom=prod(sqrt(2.*pi.*N.*z));
                    expo=sum((yi-z).^2./(z));
                    out(i,k,p,o,l)=sqrt(2.*pi.*N)/denom.*exp(-N/2.*expo);
                    %fix a slight numerical bug that can occur
                    if denom==0
                        %in this case, the exponential is zero, but dividing by zero gives
                        %a NaN, so you have make sure it is correctly recorded as a zero.
                        out(i,k,p,o,l)=0;
                    end
                end
            end
        end
    end
end
end




