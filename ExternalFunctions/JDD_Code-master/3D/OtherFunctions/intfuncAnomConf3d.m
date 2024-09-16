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

function[out]=intfuncAnomConf3d(x1,x2,x3,x4,N,yi,ri,dr,tau)
out=zeros(length(x1),length(x2),length(x3),length(x4));
for i=1:length(x2)
    Dalpha=x2(i);
    for j=1:length(x3)
        alpha=x3(j);
        if alpha < 0.5 && alpha>0.4
            min=-700; %limits on inverse laplace transform
        elseif alpha<= 0.4 && alpha>0.3
            min=-900;
        elseif alpha<=0.3
            min=-1100;
        else
            min=-700;
        end
        
        %calculate predicted probabilities
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        intsize=200;
        numint=ceil(abs(min+100)/intsize);
        z=zeros(size(ri));
        for k=1:numint
            z=z+2*dr*ri/(Dalpha).*abs(integral(fun,min+(k-1)*intsize-1i*1e-6,min+k*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        z1=z+dr*ri/(Dalpha).*abs(integral(fun,-100-1i*1e-6,100-1i*1e-6,...
            'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        for k = 1:length(x1)
            DConf = x1(k);
            z2=(4*pi*ri.^2*dr) - ((ri.^2)*dr/(2*sqrt(pi)*(DConf*tau)^(3/2))).*exp(-ri.^2/(4*DConf*tau));
            parfor l=1:length(x4)
                [i,j,k,l]
                fd=x4(l);
                z=fd*z1+(1-fd)*z2;
                %calculate P(JDD|model,parameters)
                denom=prod(sqrt(2.*pi.*N.*z));
                expo=sum((yi-z).^2./(z));
                out(i,j,k,l)=sqrt(2.*pi.*N)/denom.*exp(-N/2.*expo);
                %fix a slight numerical bug that can occur
                if denom==0
                    %in this case, the exponential is zero, but dividing by zero gives
                    %a NaN, so you have make sure it is correctly recorded as a zero.
                    out(i,j,k,l)=0;
                end
            end
        end
    end
end
end




