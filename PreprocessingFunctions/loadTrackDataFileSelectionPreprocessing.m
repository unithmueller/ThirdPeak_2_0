function TrackData = loadTrackDataFileSelectionPreprocessing(ImportSettingsStruct)
    %Function to grab a file depending on the settings set in the import
    %settings. Will allow for tracks of filtered localisation to be
    %visualized in the preprocessing window and used for further
    %calculations
    %Check that we have settings set
    %tracklength = minTrackLength;
    if isa(ImportSettingsStruct, "struct")
        %get the type of import
        imptype = ImportSettingsStruct.type;
        switch imptype
            case "Swift-AllData CSV"
                %% SwiftAllData
                % Display uigetfile dialog
                filterspec = {'*.csv'};
                [f, p] = uigetfile(filterspec, 'MultiSelect','on');
                % Make sure user didn't cancel uigetfile dialog
                if ~ischar(p)
                    error("No file selected")
                end
                if isa(f,"char") %only one file selected, convert to cell 
                    b = f;
                    f = cell(1,1);
                    f{1,1} = b;
                end
                loadfile = cell(size(f,2),2);
                for i = 1:size(f,2)
                    loadfile{i,1} = ImportSwiftAsTrc(fullfile(p,f{1,i}),1);
                    loadfile{i,2} = f{1,i};
                end
                %make sure that traces are actually loaded
                if size(loadfile{1,1},2) < 5
                    error("Not a trace file!")
                else
                    TrackData = loadfile;
                end
            case "Swift-MinData CSV" 
                %% SwiftMinData
                % Display uigetfile dialog
                filterspec = {'*.csv'};
                [f, p] = uigetfile(filterspec, 'MultiSelect','on');
                % Make sure user didn't cancel uigetfile dialog
                if ~ischar(p)
                    error("No file selected")
                end
                loadfile = cell(size(f,2),2);
                for i = 1:size(f,2)
                    loadfile(i,1) = ImportSwiftAsTrc(fullfile(p,f{1,i}),0);
                    loadfile(i,2) = f;
                end
                %make sure that traces are actually loaded
                if size(loadfile{1,1},2) < 5
                    error("Not a trace file!")
                else
                    TrackData = loadfile;
                end
            case "Custom Trackdata CSV"
                %% Custom CSV
                % Display uigetfile dialog
                filterspec = {'*.csv'};
                [f, p] = uigetfile(filterspec, 'MultiSelect','on');
                % Make sure user didn't cancel uigetfile dialog
                if ~ischar(p)
                    error("No file selected")
                end
                if isa(f,"char") %only one file selected, convert to cell 
                    b = f;
                    f = cell(1,1);
                    f{1,1} = b;
                end
                loadfile = cell(size(f,2),2);
                for i = 1:size(f,2)
                    loadfile{i,1} = loadCustomCSV(fullfile(p,f{1,i}), ImportSettingsStruct);
                    loadfile{i,2} = f{1,i};
                end
                %make sure that traces are actually loaded
                if size(loadfile{1,1},2) < 5
                    error("Not a trace file!")
                else
                    TrackData = loadfile;
                end
            case "Custom Trackdata .mat"
                %% Custom mat
                 % Display uigetfile dialog
                filterspec = {'*.mat'};
                [f, p] = uigetfile(filterspec, 'MultiSelect','on');
                % Make sure user didn't cancel uigetfile dialog
                if ~ischar(p)
                    error("No file selected")
                end
                if isa(f,"char") %only one file selected, convert to cell 
                    b = f;
                    f = cell(1,1);
                    f{1,1} = b;
                end
                loadfile = cell(size(f,2),2);
                for i = 1:size(f,2)
                    loadfile{i,1} = loadCustomMat(fullfile(p,f{1,i}), ImportSettingsStruct);
                    loadfile{i,2} = f{1,i};
                end
                %make sure that traces are actually loaded
                if size(loadfile{1,1},2) < 5
                    error("Not a trace file!")
                else
                    TrackData = loadfile;
                end
            case "Filtered data .mat"
                %% Filtered data mat
                % Display uigetfile dialog
                filterspec = {'*.mat'};
                [f, p] = uigetfile(filterspec, 'MultiSelect','on');
                % Make sure user didn't cancel uigetfile dialog
                if ~ischar(p)
                    error("No file selected")
                end
                if isa(f,"char") %only one file selected, convert to cell 
                    b = f;
                    f = cell(1,1);
                    f{1,1} = b;
                end
                loadfile = cell(size(f,2),2);
                for i = 1:size(f,2)
                    loadfile{i,1} = load(fullfile(p,f{1,i}));
                    loadfile{i,2} = f{1,i};
                end
                loadfile = loadfile{1,1}.datatosave;
                TrackData = loadfile;
        end
    end           
end