%1D integration for Bayesian Classifier
%Rebecca Menssen
%This version of the code: 3/20/19

%%%%%%%%%%INPUTS%%%%%%%%%%
%dbeta--defines the bounds on integration. Have been using 2 times the
%       standard deviation of the bootstrapped parameters.
%beta--median or mean parameter values. Center of integration window.
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%prob--vector giving normalized probabilities in the order pure diffusion,
%directed diffusion, anomalous diffusion
%value--maximal probability
%method--index giving the method of the maximal probability
% 
% 
% nameArray = ["FreeD"; "DirV";"DirDv";"AnoDalpha";"Anoalpha";"ConfD";"ConfRad";"FFD1";"FFD2";"FFfdDD";"FreeDirFreeD";"FreeDirV";"FreeDirVD";"FreeDirfdDV";
%     "FreeAnomFreeD";"FreeAnomANomD";"FreeAnomalphada";"FreeAnomfdDA"; "FreeConfFD";"FreeConfConD";"FreeConfConR";"FreeConfConfd";"DirDirD1";"DirDirD2";
%     "DirDirV1"; "DirDirV2";"DirDirfd";"DirAnomD1";"DirAnomV1";"DirAnomD2";"DirAnomAlph2";"DirAnomfd"; "DirConfD1"; "DirConfV1"; "DirConfD2"; "DirConfR2";
%     "DirConffd"; "AnomAnomD1"; "AnomAnomD2"; "AnomAnomA1"; "AnomAnomA2"; "AnomAnomfd"; "AnomConfD1"; "AnomConfA1"; "AnomConfD2"; "AnomConfR2"; "AnomConffd";
%     "ConfConfD1"; "ConfConfR1";"ConfConfD2";"ConfConfR2";"ConfConffd"];



function [prob,value,method]=ExtendedIntegration2DwithComboModels(modelsToTest,dbeta,beta,N,yi,ri,dr,tau)
%Vector for storing probabilities
prob=zeros(size(modelsToTest,2),2);

for i = 1:size(modelsToTest,1)
    switch modelsToTest(i)
        case 1%done
            %DIFFUSION INTEGRATION
            %define ranges on integration.
            minD=beta(1)-dbeta(1);
            maxD=beta(1)+dbeta(1);
            
            %Check against negative values. They tend to really mess things up and
            %should be ignored.
            if minD <=0
                minD=1e-10;
            end
            if maxD <=0
                maxD=1e-9;
            end
            
            %Find length of interval, find probabilities at different values of
            %parameters, and then do the integration.
            lengthD = maxD-minD;
            x1 = linspace(minD,maxD,500); 
            out=intfuncD2D(x1,N,yi,ri,dr,tau);
            prob(i,1) = 1;
            prob(i,2)=1/lengthD*trapz(x1,out);
        case 2%done
            %DIRECTED MOTION INTEGRATION
            %define ranges on integration.
            minV=beta(2)-dbeta(2);
            maxV=beta(2)+dbeta(2);
            minDv=beta(3)-dbeta(3);
            maxDv=beta(3)+dbeta(3);

            %Checks against negative values.
            if minV <= 0
                minV=0.01;
            end
            if maxV <= 0
                maxV=0.02;
            end
            if minDv < 0
                minDv=1e-10;
            end
            if maxDv < 0
                maxDv=1e-9;
            end

            %Find length of interval, find probabilities at different values of
            %parameters, and then do the integration.
            lengthV=maxV-minV;
            lengthDv=maxDv-minDv;
            x1 = linspace(minV,maxV,100);
            x2 = linspace(minDv,maxDv,100);
            [X1,X2] = meshgrid(x1,x2);
            out=intfuncV2D(X1,X2,N,yi,ri,dr,tau);
            prob(i,1) = 2;
            prob(i,2)=1/lengthV*1/lengthDv*trapz(x2,trapz(x1,out,1));
        case 3%done
            %ANOMALOUS DIFFUSION INTEGRATION

            minDalpha=beta(4)-dbeta(4);
            maxDalpha=beta(4)+dbeta(4);
            minalpha=beta(5)-dbeta(5);
            maxalpha=beta(5)+dbeta(5);

            %checks against negative values
            if minDalpha <=0
                minDalpha=1e-10;
            end
            if maxDalpha <=0
                maxDalpha=1e-9;
            end
            if minalpha <= 0
                minalpha=0.01;
            end
            if maxalpha <= 0
                maxalpha=0.02;
            end

            %now the case where alpha>1
            prob(i,1) = 3;
            if maxalpha>1 && minalpha>1
                prob(i,2)=0;
            elseif maxalpha>1 %just upper bound is above.
                maxalpha=1;
                lengthDalpha=maxDalpha-minDalpha;
                lengthalpha=maxalpha-minalpha;
                x1 = linspace(minDalpha,maxDalpha,100);
                x2 = linspace(minalpha,maxalpha,100);
                [X1,X2] = meshgrid(x1,x2);
                out=intfuncA2D(X1,X2,N,yi,ri,dr,tau);
                prob(i,2)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));
            else %both below 1.
                lengthDalpha=maxDalpha-minDalpha;
                lengthalpha=maxalpha-minalpha;
                x1 = linspace(minDalpha,maxDalpha,100);
                x2 = linspace(minalpha,maxalpha,100);
                [X1,X2] = meshgrid(x1,x2);
                out=intfuncA2D(X1,X2,N,yi,ri,dr,tau);
                prob(i,2)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));
            end
        case 4%done
            %confined
            minDalpha=beta(6)-dbeta(6);
            maxDalpha=beta(6)+dbeta(6);
            %checks against negative values
            if minDalpha <=0
                minDalpha=1e-10;
            end
            if maxDalpha <=0
                maxDalpha=1e-9;
            end
            prob(i,1) = 4;
            lengthDalpha=maxDalpha-minDalpha;
            x1 = linspace(minDalpha,maxDalpha,100);
            [out]=intfuncConf2d(x1,N,yi,ri,dr,tau);
            prob(i,2)=1/lengthDalpha*trapz(x1,out,1);
        case 5%done
            %%double diffusion
            minD1=beta(8)-dbeta(8); maxD1=beta(8)+dbeta(8);
            minD2=beta(9)-dbeta(9); maxD2=beta(9)+dbeta(9);
            minfd=beta(10)-dbeta(10); maxfd=beta(10)+dbeta(10);
            if minD1 <=0, minD1=1e-10; end
            if maxD1 <=0, maxD1=1e-9; end
            if minD2 <=0, minD2=1e-10; end
            if maxD2 <=0, maxD2=1e-9; end
            if minfd <=0, minfd=1e-10; end
            if maxfd <=0, maxfd=1e-9; end
            if maxfd > 1, maxfd=1; end
            lengthD1=maxD1-minD1; lengthD2=maxD2-minD2; lengthfd=maxfd-minfd;
            x1 = linspace(minD1,maxD1,50); x2 = linspace(minD2,maxD2,50); x3 = linspace(minfd,maxfd,24);
            [out]=intfuncDD2d(x1,x2,x3,N,yi,ri,dr,tau);
            prob(i,1) = 5;
            prob(i,2)=1/lengthD1*1/lengthD2*1/lengthfd*trapz(x3,(trapz(x2,trapz(x1,out,1),2)));
        case 6%done
            %%diffusion directed combo motion
            minDdv=beta(11)-dbeta(11); maxDdv=beta(11)+dbeta(11);
            minVdv=beta(12)-dbeta(12); maxVdv=beta(12)+dbeta(12);
            minDvdv=beta(13)-dbeta(13); maxDvdv=beta(13)+dbeta(13);
            minfd=beta(14)-dbeta(14); maxfd=beta(14)+dbeta(14);
            if minDdv <=0, minDdv=1e-10; end
            if maxDdv <=0, maxDdv=1e-9; end
            if minVdv <=0, minVdv=1e-10; end
            if maxVdv <=0, maxVdv=1e-9; end
            if minDvdv <=0, minDvdv=1e-10; end
            if maxDvdv <=0, maxDvdv=1e-9; end
            if minfd <=0, minfd=1e-10; end
            if maxfd <=0, maxfd=1e-9; end
            if maxfd > 1, maxfd=1; end
            lengthDdv=maxDdv-minDdv; lengthVdv=maxVdv-minVdv; lengthDvdv=maxDvdv-minDvdv; lengthfd=maxfd-minfd;
            x1 = linspace(minDdv,maxDdv,12); x2 = linspace(minVdv,maxVdv,12);
            x3 = linspace(minDvdv,maxDvdv,12); x4 = linspace(minfd,maxfd,12);
            out=intfuncDV2d(x1,x2,x3,x4,N,yi,ri,dr,tau);
            prob(i,1) = 6;
            prob(i,2)=1/lengthDdv*1/lengthVdv*1/lengthDvdv*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));
        case 7
            %%diffusion anomalous combo motion
            minDda=beta(15)-dbeta(15); maxDda=beta(15)+dbeta(15);
            minDalphada=beta(16)-dbeta(16); maxDalphada=beta(16)+dbeta(16);
            minalphada=beta(17)-dbeta(17); maxalphada=beta(17)+dbeta(17);
            minfd=beta(18)-dbeta(18); maxfd=beta(18)+dbeta(18);
            if minDda <=0, minDda=1e-14; end
            if maxDda <=0, maxDda=1e-13; end
            if minDalphada <=0, minDalphada=1e-14; end
            if maxDalphada <=0, maxDalphada=1e-13; end
            if minalphada <=0, minalphada=1e-14; end
            if maxalphada <=0, maxalphada=1e-13; end
            if minfd <=0, minfd=1e-14; end
            if maxfd <=0, maxfd=1e-13; end
            if maxfd > 1, maxfd=1; end
            %now cases where alpha is greater than 1
            prob(i,1) = 7;
            if maxalphada && minalphada > 1
                prob(i,2)=0;
            elseif maxalphada>1
                maxalphada=1;
                lengthDda=maxDda-minDda; lengthDalphada=maxDalphada-minDalphada;
                lengthalphada=maxalphada-minalphada; lengthfd=maxfd-minfd;
                x1 = linspace(minDda,maxDda,24); x2 = linspace(minDalphada,maxDalphada,24);
                x3 = linspace(minalphada,maxalphada,24); x4 = linspace(minfd,maxfd,24);
                [out]=intfuncDA2d(x1,x2,x3,x4,N,yi,ri,dr,tau);
                prob(i,2)=1/lengthDda*1/lengthDalphada*1/lengthalphada*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));
            else
                lengthDda=maxDda-minDda; lengthDalphada=maxDalphada-minDalphada;
                lengthalphada=maxalphada-minalphada; lengthfd=maxfd-minfd;
                x1 = linspace(minDda,maxDda,12); x2 = linspace(minDalphada,maxDalphada,12);
                x3 = linspace(minalphada,maxalphada,12); x4 = linspace(minfd,maxfd,12);
                [out]=intfuncDA2d(x1,x2,x3,x4,N,yi,ri,dr,tau);
                prob(i,2)=1/lengthDda*1/lengthDalphada*1/lengthalphada*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));
            end
        case 8%done
            %freeconfined
            %define ranges on integration.
            minD=beta(19)-dbeta(19);
            maxD=beta(19)+dbeta(19);
            minDConf=beta(20)-dbeta(20);
            maxDConf=beta(20)+dbeta(20);
            minfd=beta(22)-dbeta(22); maxfd=beta(22)+dbeta(22);
            %Check against negative values. They tend to really mess things up and
            %should be ignored.
            if minD <=0
                minD=1e-10;
            end
            if maxD <=0
                maxD=1e-9;
            end
            if minDConf <=0
                minDConf=1e-10;
            end
            if maxDConf <=0
                maxDConf=1e-9;
            end
            if minfd <=0, minfd=1e-14; end
            if maxfd <=0, maxfd=1e-13; end
            if maxfd > 1, maxfd=1; end
            
            %Find length of interval, find probabilities at different values of
            %parameters, and then do the integration.
            lengthD = maxD-minD;
            lengthDConf=maxDConf-minDConf;
            lengthfd=maxfd-minfd;

            x1 = linspace(minD,maxD,50);
            x2 = linspace(minDConf,maxDConf,50);
            x3 = linspace(minfd,maxfd,50);
            
            %checks against negative values

            prob(i,1) = 8;
            [out]=intfuncDConf2d(x1,x2,x3,N,yi,ri,dr,tau);
            prob(i,2)=1/lengthD*1/lengthDConf*1/lengthfd*trapz(x3,(trapz(x2,trapz(x1,out,1),2)));
        case 9%done
            %directdirect
            minVdv1=beta(25)-dbeta(25); maxVdv1=beta(25)+dbeta(25);
            minDvdv1=beta(23)-dbeta(23); maxDvdv1=beta(23)+dbeta(23);
            minVdv2=beta(26)-dbeta(26); maxVdv2=beta(26)+dbeta(26);
            minDvdv2=beta(24)-dbeta(24); maxDvdv2=beta(24)+dbeta(24);
            minfd=beta(27)-dbeta(27); maxfd=beta(27)+dbeta(27);
            if minVdv1 <=0, minVdv1=1e-10; end
            if maxVdv1 <=0, maxVdv1=1e-9; end
            if minDvdv1 <=0, minDvdv1=1e-10; end
            if maxDvdv1 <=0, maxDvdv1=1e-9; end
            if minVdv2 <=0, minVdv2=1e-10; end
            if maxVdv2 <=0, maxVdv2=1e-9; end
            if minDvdv2 <=0, minDvdv2=1e-10; end
            if maxDvdv2 <=0, maxDvdv2=1e-9; end
            if minfd <=0, minfd=1e-10; end
            if maxfd <=0, maxfd=1e-9; end
            if maxfd > 1, maxfd=1; end

            lengthDdv1=maxDvdv1-minDvdv1; lengthVdv1=maxVdv1-minVdv1; lengthDvdv2=maxDvdv2-minDvdv2; lengthVdv2=maxVdv2-minVdv2; lengthfd=maxfd-minfd;
            x1 = linspace(minDvdv1,maxDvdv1,12); x2 = linspace(minVdv1,maxVdv1,12);
            x3 = linspace(minDvdv2,maxDvdv2,12); x4 = linspace(minVdv2,maxVdv2,12);
            x5 = linspace(minfd,maxfd,12);
            [out]=intfuncDirDir2d(x1,x2,x3,x4,x5,N,yi,ri,dr,tau);
            prob(i,1) = 9;
            prob(i,2)=1/lengthDdv1*1/lengthVdv1*1/lengthDvdv2*1/lengthVdv2*1/lengthfd*trapz(x5,trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3),4));
        case 10
            %directanomal%done
            minVdv1=beta(29)-dbeta(29); maxVdv1=beta(29)+dbeta(29);
            minDvdv1=beta(28)-dbeta(28); maxDvdv1=beta(28)+dbeta(28);
            minDA=beta(30)-dbeta(30); maxDA=beta(30)+dbeta(30);
            minAlph=beta(31)-dbeta(31); maxAlph=beta(31)+dbeta(31);
            minfd=beta(32)-dbeta(32); maxfd=beta(32)+dbeta(32);
            if minVdv1 <=0, minVdv1=1e-10; end
            if maxVdv1 <=0, maxVdv1=1e-9; end
            if minDvdv1 <=0, minDvdv1=1e-10; end
            if maxDvdv1 <=0, maxDvdv1=1e-9; end
            if minDA <=0, minDA=1e-10; end
            if maxDA <=0, maxDA=1e-9; end
            if minAlph <=0, minAlph=1e-10; end
            if maxAlph <=0, maxAlph=1e-9; end
            if minfd <=0, minfd=1e-10; end
            if maxfd <=0, maxfd=1e-9; end
            if maxfd > 1, maxfd=1; end

            lengthDdv1=maxDvdv1-minDvdv1; lengthVdv1=maxVdv1-minVdv1; lengthDA=maxDA-minDA; lengthAlph=maxAlph-minAlph; lengthfd=maxfd-minfd;
            x1 = linspace(minDvdv1,maxDvdv1,12); x2 = linspace(minVdv1,maxVdv1,12);
            x3 = linspace(minAlph,maxAlph,12); x4 = linspace(minDA,maxDA,12);
            x5 = linspace(minfd,maxfd,12);
            [out]=intfuncVA2D(x2,x1,x4,x3,x5,N,yi,ri,dr,tau);
            prob(i,1) = 10;
            prob(i,2)=1/lengthDdv1*1/lengthVdv1*1/lengthDA*1/lengthAlph*1/lengthfd*trapz(x5,trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3),4));
        case 11
            %directconfined
            minVdv1=beta(34)-dbeta(34); maxVdv1=beta(34)+dbeta(34);
            minDvdv1=beta(33)-dbeta(33); maxDvdv1=beta(33)+dbeta(33);
            minDConf=beta(35)-dbeta(35); maxDConf=beta(35)+dbeta(35);
            minfd=beta(37)-dbeta(37); maxfd=beta(37)+dbeta(37);
            if minVdv1 <=0, minVdv1=1e-10; end
            if maxVdv1 <=0, maxVdv1=1e-9; end
            if minDvdv1 <=0, minDvdv1=1e-10; end
            if maxDvdv1 <=0, maxDvdv1=1e-9; end
            if minDConf <=0, minDConf=1e-10; end
            if maxDConf <=0, maxDConf=1e-9; end
            if minfd <=0, minfd=1e-10; end
            if maxfd <=0, maxfd=1e-9; end
            if maxfd > 1, maxfd=1; end

            lengthV=maxVdv1-minVdv1; lengthDV=maxDvdv1-minDvdv1;
            lengthDConf=maxDConf-minDConf; lengthfd=maxfd-minfd;
            prob(i,1) = 11;
            x1 = linspace(minVdv1,maxVdv1,12); x2 = linspace(minDvdv1,maxDvdv1,12);
            x3 = linspace(minDConf,maxDConf,12); x4 = linspace(minfd,maxfd,12);
            [out]=intfuncDirConf2d(x1,x2,x3,x4,N,yi,ri,dr,tau);
            prob(i,2)=1/lengthV*1/lengthDV*1/lengthDConf*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));

        case 12
            %anomal anomal
            minDalpha1=beta(38)-dbeta(38);
            maxDalpha1=beta(38)+dbeta(38);
            minalpha1=beta(40)-dbeta(40);
            maxalpha1=beta(40)+dbeta(40);
            minDalpha2=beta(39)-dbeta(39);
            maxDalpha2=beta(39)+dbeta(39);
            minalpha2=beta(41)-dbeta(41);
            maxalpha2=beta(41)+dbeta(41);
            minfd=beta(42)-dbeta(42); maxfd=beta(42)+dbeta(42);
            %checks against negative values
            if minDalpha1 <=0, minDalpha1=1e-10; end
            if maxDalpha1 <=0,maxDalpha1=1e-9; end
            if minalpha1 <= 0, minalpha1=0.01; end
            if maxalpha1 <= 0, maxalpha1=0.02; end
            if minDalpha2 <=0, minDalpha2=1e-10; end
            if maxDalpha2 <=0,maxDalpha2=1e-9; end
            if minalpha2 <= 0, minalpha2=0.01; end
            if maxalpha2 <= 0, maxalpha2=0.02; end
            if minfd <=0, minfd=1e-10; end
            if maxfd <=0, maxfd=1e-9; end
            if maxfd > 1, maxfd=1; end

            %now the case where alpha>1
            prob(i,1) = 12;
            if maxalpha1>1 && minalpha1>1 && maxalpha2>1 && minalpha2>1
                prob(i,2)=0;
            elseif maxalpha1>1 %just upper bound is above.
                maxalpha1=1;
                if maxalpha2>1
                    maxalpha2 = 1;
                end
                lengthDalpha1=maxDalpha1-minDalpha1;
                lengthalpha1=maxalpha1-minalpha1;
                lengthDalpha2=maxDalpha2-minDalpha2;
                lengthalpha2=maxalpha2-minalpha2;
                lengthfd=maxfd-minfd;

                x1 = linspace(minDalpha1,maxDalpha1,12);
                x2 = linspace(minalpha1,maxalpha1,12);
                x3 = linspace(minDalpha2,maxDalpha2,12);
                x4 = linspace(minalpha2,maxalpha2,12);
                x5 = linspace(minfd,maxfd,12);

                [out]=intfuncAnoAno2d(x1,x2,x3,x4, x5, N,yi,ri,dr,tau);
                prob(i,2)=1/lengthDalpha1*1/lengthalpha1*1/lengthDalpha2*1/lengthalpha2*1/lengthfd*trapz(x5,trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3),4));

            else %both below 1.
                if maxalpha2>1
                    maxalpha2 = 1;
                end
                lengthDalpha1=maxDalpha1-minDalpha1;
                lengthalpha1=maxalpha1-minalpha1;
                lengthDalpha2=maxDalpha2-minDalpha2;
                lengthalpha2=maxalpha2-minalpha2;
                lengthfd=maxfd-minfd;

                x1 = linspace(minDalpha1,maxDalpha1,12);
                x2 = linspace(minalpha1,maxalpha1,12);
                x3 = linspace(minDalpha2,maxDalpha2,12);
                x4 = linspace(minalpha2,maxalpha2,12);
                x5 = linspace(minfd,maxfd,12);

                [out]=intfuncAnoAno2d(x1,x2,x3,x4, x5, N,yi,ri,dr,tau);
                prob(i,2)=1/lengthDalpha1*1/lengthalpha1*1/lengthDalpha2*1/lengthalpha2*1/lengthfd*trapz(x5,trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3),4));
            end

        case 13
            %anomal confined
            minDConf=beta(45)-dbeta(45); maxDConf=beta(45)+dbeta(45);
            minDalphada=beta(43)-dbeta(43); maxDalphada=beta(43)+dbeta(43);
            minalphada=beta(44)-dbeta(44); maxalphada=beta(44)+dbeta(44);
            minfd=beta(47)-dbeta(47); maxfd=beta(47)+dbeta(47);
            if minDConf <=0, minDConf=1e-14; end
            if maxDConf <=0, maxDConf=1e-13; end
            if minDalphada <=0, minDalphada=1e-14; end
            if maxDalphada <=0, maxDalphada=1e-13; end
            if minalphada <=0, minalphada=1e-14; end
            if maxalphada <=0, maxalphada=1e-13; end
            if minfd <=0, minfd=1e-14; end
            if maxfd <=0, maxfd=1e-13; end
            if maxfd > 1, maxfd=1; end
            %now cases where alpha is greater than 1
            prob(i,1) = 13;
            if maxalphada && minalphada > 1
                prob(i,2)=0;
            elseif maxalphada>1
                maxalphada=1;
                lengthDConf=maxDConf-minDConf; lengthDalphada=maxDalphada-minDalphada;
                lengthalphada=maxalphada-minalphada; lengthfd=maxfd-minfd;
                x1 = linspace(minDConf,maxDConf,24); x2 = linspace(minDalphada,maxDalphada,24);
                x3 = linspace(minalphada,maxalphada,24); x4 = linspace(minfd,maxfd,24);
                [out]=intfuncAnomConf2d(x1,x2,x3,x4,N,yi,ri,dr,tau);
                prob(i,2)=1/lengthDConf*1/lengthDalphada*1/lengthalphada*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));
            else
                lengthDConf=maxDConf-minDConf; lengthDalphada=maxDalphada-minDalphada;
                lengthalphada=maxalphada-minalphada; lengthfd=maxfd-minfd;
                x1 = linspace(minDConf,maxDConf,12); x2 = linspace(minDalphada,maxDalphada,12);
                x3 = linspace(minalphada,maxalphada,12); x4 = linspace(minfd,maxfd,12);
                [out]=intfuncAnomConf2d(x1,x2,x3,x4,N,yi,ri,dr,tau);
                prob(i,2)=1/lengthDConf*1/lengthDalphada*1/lengthalphada*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));
            end
        case 14
            %confined confined
            minDConf1=beta(6)-dbeta(6); maxDConf1=beta(6)+dbeta(6);
            minDConf2=beta(6)-dbeta(6); maxDConf2=beta(6)+dbeta(6);
            minfd=beta(47)-dbeta(47); maxfd=beta(47)+dbeta(47);
            if minDConf1 <=0, minDConf1=1e-14; end
            if maxDConf1 <=0, maxDConf1=1e-13; end
            if minDConf2 <=0, minDConf2=1e-14; end
            if maxDConf2 <=0, maxDConf2=1e-13; end
            if minfd <=0, minfd=1e-14; end
            if maxfd <=0, maxfd=1e-13; end
            if maxfd > 1, maxfd=1; end
            %checks against negative values
            prob(i,1) = 14;
            lengthDConf1=maxDConf1-minDConf1;
            lengthDConf2=maxDConf2-minDConf2;
            lengthfd=maxfd-minfd;
            x1 = linspace(minDConf1,maxDConf1,50);
            x2 = linspace(minDConf2,maxDConf2,50);
            x3 = linspace(minfd,maxfd,50);
            [out]=intfuncConfConf2d(x1,x2,x3,N,yi,ri,dr,tau);
            prob(i,2)=1/lengthDConf1*1/lengthDConf2*1/lengthfd*trapz(x3,(trapz(x2,trapz(x1,out,1),2)));
    end
end


%Normalize Probabilities (if wanted, not necessary really). 
prob(:,2)=prob(:,2)./sum(prob(:,2));

%Select the best method
value = max(prob(:,2));
method = prob(prob(:,2)==value,1);
end

