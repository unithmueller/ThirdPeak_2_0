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

function[out]=intfuncDirAnom(x1,x2,x3,x4,x5,N,yi,ri,dr,tau)
out=zeros(length(x1),length(x2),length(x3),length(x4),length(x5));
for i = 1:length(x2)
    V1 = x2(i);
    for p=1:length(x1)
        Dv1=x1(p);
        z1 = dr/((4*pi*Dv1*tau)^(1/2)).*exp(-(ri.^2+V1^2*tau^2)/(4*Dv1*tau)+ri*V1/(2*Dv1))+...
            dr/((4*pi*Dv1*tau)^(1/2)).*exp(-(ri.^2+V1^2*tau^2)/(4*Dv1*tau)-ri*V1/(2*Dv1));
        for j=1:length(x4)
            Dalpha=x4(j);
            for k=1:length(x3)
                alpha=x3(k);
                if alpha < 0.5
                    min=-300^(.5/alpha); %limits on inverse laplace transform
                else
                    min=-500;
                end
                %calculate predicted probabilities
                fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*...
                    exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
                z2=dr/(Dalpha)^(1/2).*abs(integral(fun,min,-1*min,...
                    'ArrayValued',true,'AbsTol',2e-5));
                parfor l=1:length(x5)
                    [i,p,j,k,l]
                    fd=x5(l);
                    z=fd*z1+(1-fd)*z2;
                    %calculate P(JDD|model,parameters)
                    denom=prod(sqrt(2.*pi.*N.*z));
                    expo=sum((yi-z).^2./(z));
                    out(i,p,j,k,l)=sqrt(2.*pi.*N)/denom.*exp(-N/2.*expo);
                    %fix a slight numerical bug that can occur
                    if denom==0
                        %in this case, the exponential is zero, but dividing by zero gives
                        %a NaN, so you have make sure it is correctly recorded as a zero.
                        out(i,p,j,k,l)=0;
                    end
                end
            end
        end
    end
end
end




