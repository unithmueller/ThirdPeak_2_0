function [strucNames, nameArray, indArray] = getSaveStrucNamesJDD(diffModelUsed)
%Returns names of the save structure properties in dependency of the models
%used
output = {};

for i = 1:size(diffModelUsed,1)
    switch diffModelUsed(i)
        case 1
            output{i,1} = "FreeD";
        case 2
            output{i,1}  = "DirV"; output{i,2}  = "DirDv";
        case 3
            output{i,1}  = "AnoDalpha"; output{i,2} = "Anoalpha";
        case 4
            output{i,1}  = "ConfD"; output{i,2} = "ConfRad";
        case 5
            output{i,1}  = "FFD1"; output{i,2} = "FFD2"; output{i,3} = "FFfdDD";
        case 6
            output{i,1}  = "FreeDirFreeD"; output{i,2} = "FreeDirV"; output{i,3} = "FreeDirVD"; output{i,4} = "FreeDirfdDV";
        case 7
            output{i,1}  = "FreeAnomFreeD"; output{i,2} = "FreeAnomAnomD"; output{i,3} = "FreeAnomalphada"; output{i,4} = "FreeAnomfdDA";
        case 8
            output{i,1}  = "FreeConfFD"; output{i,2} = "FreeConfConD"; output{i,3} = "FreeConfConR"; output{i,4} = "FreeConfConfd";
        case 9
            output{i,1}  = "DirDirD1"; output{i,2} = "DirDirD2"; output{i,3} = "DirDirV1"; output{i,4} = "DirDirV2"; output{i,5} = "DirDirfd";
        case 10
            output{i,1}  = "DirAnomD1"; output{i,2} = "DirAnomV1"; output{i,3} = "DirAnomD2"; output{i,4} = "DirAnomAlph2"; output{i,5} = "DirAnomfd";
        case 11
            output{i,1}  = "DirConfD1"; output{i,2} = "DirConfV1"; output{i,3} = "DirConfD2"; output{i,4} = "DirConfR2"; output{i,5} = "DirConffd";
        case 12
            output{i,1}  = "AnomAnomD1"; output{i,2} = "AnomAnomD2"; output{i,3} = "AnomAnomA1"; output{i,4} = "AnomAnomA2"; output{i,5} = "AnomAnomfd";
        case 13
            output{i,1}  = "AnomConfD1"; output{i,2} = "AnomConfA1"; output{i,3} = "AnomConfD2"; output{i,4} = "AnomConfR2"; output{i,5} = "AnomConffd";
        case 14
            output{i,1}  = "ConfConfD1"; output{i,2} = "ConfConfR1"; output{i,3} = "ConfConfD2"; output{i,4} = "ConfConfR2"; output{i,5} = "ConfConffd";
    end
end

strucNames = [output{:}];

nameArray = ["FreeD"; "DirV";"DirDv";"AnoDalpha";"Anoalpha";"ConfD";"ConfRad";"FFD1";"FFD2";"FFfdDD";"FreeDirFreeD";"FreeDirV";"FreeDirVD";"FreeDirfdDV";
    "FreeAnomFreeD";"FreeAnomAnomD";"FreeAnomalphada";"FreeAnomfdDA"; "FreeConfFD";"FreeConfConD";"FreeConfConR";"FreeConfConfd";"DirDirD1";"DirDirD2";
    "DirDirV1"; "DirDirV2";"DirDirfd";"DirAnomD1";"DirAnomV1";"DirAnomD2";"DirAnomAlph2";"DirAnomfd"; "DirConfD1"; "DirConfV1"; "DirConfD2"; "DirConfR2";
    "DirConffd"; "AnomAnomD1"; "AnomAnomD2"; "AnomAnomA1"; "AnomAnomA2"; "AnomAnomfd"; "AnomConfD1"; "AnomConfA1"; "AnomConfD2"; "AnomConfR2"; "AnomConffd";
    "ConfConfD1"; "ConfConfR1";"ConfConfD2";"ConfConfR2";"ConfConffd"];

indArray = zeros(size(strucNames,2),1);

for i = 1:size(strucNames,2)
    indArray(i) = find(nameArray == strucNames(i));
end

end