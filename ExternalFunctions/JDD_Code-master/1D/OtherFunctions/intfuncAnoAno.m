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

function[out]=intfuncAnoAno(x1,x2,x3,x4, x5, N,yi,ri,dr,tau)
out=zeros(length(x1),length(x2),length(x3),length(x4),length(x5));
for i=1:length(x1)
    Dalpha1=x1(i);
    for k=1:length(x2)
        alpha1=x2(k);
        if alpha1 < 0.5
            min1=-300^(.5/alpha1); %limits on inverse laplace transform
        else
            min1=-500;
        end
        %calculate predicted probabilities
        fun1=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha1/2-1)/(2.*pi).*...
            exp(-ri./(sqrt(Dalpha1)).*((1i*p)^(alpha1/2)));
        z1=dr/(Dalpha1)^(1/2).*abs(integral(fun1,min1,-1*min1,...
            'ArrayValued',true,'AbsTol',2e-5));
        for p=1:length(x3)
            Dalpha2=x3(p);
            for o=1:length(x4)
                alpha2=x4(o);
                if alpha2 < 0.5
                    min2=-300^(.5/alpha); %limits on inverse laplace transform
                else
                    min2=-500;
                end
                %calculate predicted probabilities
                fun2=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha2/2-1)/(2.*pi).*...
                    exp(-ri./(sqrt(Dalpha2)).*((1i*p)^(alpha2/2)));
                z2=dr/(Dalpha2)^(1/2).*abs(integral(fun2,min2,-1*min2,...
                    'ArrayValued',true,'AbsTol',2e-5));

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




