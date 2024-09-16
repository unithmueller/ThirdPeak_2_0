function pathString = SettingPythonEnvironment()
    %% Set the python environment
    %Loads the text file with the path to the python environment countaining
    %DeepSPT-packages
    
    %% 
    %Grab the file
    functionName = 'SettingPythonEnvironment.m';
    pythonPathFilePath = which(functionName);
    pythonPathFilePath = pythonPathFilePath(1:end-(size(functionName,2)));
    pythonPathEnvPath = join([pythonPathFilePath "PythonDeepSPTEnvironmentPath.txt"],"");
    pythonPathDeepSPTPath = join([pythonPathFilePath "DeepSPTLoactionPath.txt"],"");
    
    %Read the file
    pythonPathText = readlines(pythonPathEnvPath);
    pythonDeepSPTText = readlines(pythonPathDeepSPTPath);

    %put the path into a list as it make some problems with the
    %delimiter...
    splitDeepSPTPath = split(pythonDeepSPTText,["/", "\"]);
    newPath = join(splitDeepSPTPath,"/");
    
    %% Set the python environment
    try
        terminate(pyenv)
    catch
    end
    try
        %Change ExecutionMode to InProcess once everything is sorted out...
        pe = pyenv(Version=pythonPathText, ExecutionMode="OutOfProcess");
        py.os.chdir(newPath);
        pathString = py.os.getcwd();
    catch
        pathString = "";
    end
        
end

