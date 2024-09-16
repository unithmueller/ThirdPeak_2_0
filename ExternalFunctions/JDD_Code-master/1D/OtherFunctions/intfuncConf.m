%Integration function for 1D Anomalous Diffusion
%Rebecca Menssen
%Last Updated: 4/9/19

%This function creates the grid necessary to integrate across a range of
%anomalous diffusion parameters using trapz

%%%%%%%%%%INPUTS%%%%%%%%%%
%x1,x2--equally spaced array of Dalpha and alpha values to integrate over
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%out--structure of probability values for a range of Dalpha and alpha
%values

function[out]=intfuncConf(X1,N,yi,ri,dr,tau)
out=zeros(length(X1),1);
parfor i = 1:length(X1)
    %helpful to know where you are, can comment out
    DConf=X1(i);

    %calculate predicted probabilities
    z=dr/((pi*DConf*tau)^(1/2)).*exp(-ri.^2/(4*DConf*tau));

    %calculate P(JDD|model,parameters)
    denom=prod(sqrt(2*pi*N*z));
    expo=sum((yi-z).^2./(z));
    out(i)=sqrt(2*pi*N)/denom*exp(-N/2*expo);

    %fix a slight numerical bug that can occur
    if denom==0
        %in this case, the exponential is zero, but dividing by zero gives
        %a NaN, so you have make sure it is correctly recorded as a zero.
        out(i)=0;
    end
    if exp(-N/2.*expo)==0
        %in this case, the exponential is zero, but occasionally things can
        %go a bit haywire.
        out(i)=0;
    end
end
end