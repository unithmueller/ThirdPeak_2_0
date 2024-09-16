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

function[out]=intfuncDirDir3d(x1,x2,x3,x4,x5,N,yi,ri,dr,tau)
out=zeros(length(x1),length(x2),length(x3),length(x4),length(x5));
for i = 1:length(x2)
    V1 = x2(i);
    for p=1:length(x1)
        Dv1=x1(p);
        a1 = -(ri.^2+V1^2*tau^2)/(4*Dv1*tau);
        b1 = ri*V1/(2*Dv1);
        z1=dr*ri.^2*4*pi/((4*pi*Dv1*tau)^(3/2)).*exp(a1).*sinh(b1)./b1;
        for j=1:length(x4)
            V=x4(j);
            for k=1:length(x3)
                Dv=x3(k);
                a2 = -(ri.^2+V^2*tau^2)/(4*Dv*tau);
                b2 = ri*V/(2*Dv);
                z2=dr*ri.^2*4*pi/((4*pi*Dv*tau)^(3/2)).*exp(a2).*sinh(b2)./b2;
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




