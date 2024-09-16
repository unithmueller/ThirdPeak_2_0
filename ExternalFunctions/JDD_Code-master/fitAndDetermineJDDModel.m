function [param, unpackedParams, unpackedBootParams, prob, value, method, strucNames, nameArray, indArray, diffModes] = fitAndDetermineJDDModel(selectedDiffusionModes, windowType, bootstrapIterations, dimension, analyze3dAs2d, tau, dr, ri, yi, Ni, N, points, dt, calcuJumpDists,x,y,z, Nb, doModelFitting)
    %Function to manage the fitting of the different diffusion modes
    a = 1;
    %% Determine dimensionality
    if size(dimension,2) == 1
        dimVal = 1;
        x2 = [];
        x3 = [];
        if dimension == 3
            x1 = x;
        elseif dimension == 4
            x1 = y;
        elseif dimension == 5
            x1 = z;
        end
    elseif size(dimension,2) == 2
        dimVal = 2;
        x1 = x; x2 = y; x3 = [];
    elseif size(dimension,2) == 3
        dimVal = 3;
        x1 = x; x2=y; x3=z;
    end
    if analyze3dAs2d
        dimVal = 2;
        x1 = x; x2=y;
    end
    

    %% Determine the diffusion modes to fit
    diffModes = zeros(15,1);
    for i = 1:size(selectedDiffusionModes)
        diffModes(i) = selectedDiffusionModes(i).NodeData;
    end
    diffModes(diffModes(:) == 0) = [];

    %% Model fitting
    %Divide first by dimensionality, then by diffusion type
    %1 = FreeDiff   6 = Free+Directed      11 = Directed+Confined
    %2 = Directed   7 = Free+Anomal        12 = Anomal+Anomal
    %3 = Anomalous  8 = Free+Confined      13 = Anomal+Confined
    %4 = Confined   9 = Directed+Directed  14 = Confined+Confined
    %5 = Free+Free  10 = Directed+Anomal
    if dimVal == 1
        %1D analysis
        param = ExtendedModelFitting1DwithComboModels(diffModes, windowType, tau, dr, ri, yi, Ni ,N, points, dt, x1);
    elseif dimVal == 2
        %2D analysis
        param = ExtendedModelFitting2DwithComboModels(diffModes, windowType, tau, dr, ri, yi, Ni ,N, points, dt, x1, x2);
    elseif dimVal == 3
        %3D analysis
        param = ExtendedModelFitting3DwithComboModels(diffModes, windowType, tau, dr, ri, yi, Ni ,N, points, dt, x1, x2, x3);
    end

    %% Bootstrapping
    bootStructs = {};

    parfor i=1:bootstrapIterations
        X = randi(N,N,1);
        jdB=calcuJumpDists(X);
        [drB, NiB, yiB, riB] =  BinningHist(jdB, N, Nb,'no');
        if dimVal == 1
            %1D analysis
            paramB = ExtendedModelFitting1DwithComboModels(diffModes, windowType, tau, drB, riB, yiB, NiB ,N, points, dt, x1);
        elseif dimVal == 2
            %2D analysis
            paramB = ExtendedModelFitting2DwithComboModels(diffModes, windowType, tau, dr, ri, yi, Ni ,N, points, dt, x1, x2);
        elseif dimVal == 3
            %3D analysis
            paramB = ExtendedModelFitting3DwithComboModels(diffModes, windowType, tau, dr, ri, yi, Ni ,N, points, dt, x1, x2, x3);
        end
        bootStructs{i} = paramB;
    end

    %% Unpacking
    [strucNames, nameArray, indArray] = getSaveStrucNamesJDD(diffModes);
    unpackedParams = zeros(size(nameArray,1),1);
    for i = 1:size(strucNames,2)
        unpackedParams(indArray(i),1) = param.(strucNames(i));
    end
    unpackedBootParams = zeros(size(nameArray,1),bootstrapIterations);
    for i = 1:bootstrapIterations
        tmpS = bootStructs{i};
        for j = 1:size(strucNames,2)
            unpackedBootParams(indArray(j),i) = tmpS.(strucNames(j));
        end
    end
    unpackedBootParams = 2*std(unpackedBootParams,0,2);

    %% Model selection
    if doModelFitting
        if dimVal == 1
            %1D analysis
            [prob,value,method]=ExtendedIntegration1DwithComboModels(diffModes,unpackedBootParams,unpackedParams,N,yi,ri,dr,tau);
        elseif dimVal == 2
            %2D analysis
            [prob,value,method]=ExtendedIntegration2DwithComboModels(diffModes,unpackedBootParams,unpackedParams,N,yi,ri,dr,tau);
        elseif dimVal == 3
            %3D analysis
            [prob,value,method]=ExtendedIntegration3DwithComboModels(diffModes,unpackedBootParams,unpackedParams,N,yi,ri,dr,tau);
        end
    else
        prob = ones(size(diffModes,1),2); prob(:,1) = diffModes; value=1; method=1;
    end


end