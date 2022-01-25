function [ acqAlgoPprocSignals, rxData] = pre_process(sdrParams)
%%% pre_process is responsible for configuration of settings
% and pre-conditioning of the input data.
% This includes per acquisition algorithm parameters
% and common signal preparations.

%%% Initialize the settings

% Read data from configured channels 
% from current file.
rxData = read_file_data(sdrParams);


print_string('Generating C/A code.');

% Gen/Read CA/Code
caCodeTable = gen_ca_code(sdrParams.stateParams.dataPathIn, ...
                          sdrParams.sysParams.acqSatelliteList);

% Prepare signals for each algorithm
          
acqAlgoPprocSignals = cell(0);
acqAlgosList = sdrParams.sysParams.acqAlgosList;
for acqIdx = 1:length(acqAlgosList)
    
    algoName = acqAlgosList{acqIdx};    
    
    switch algoName
        case 'norm_acq_parcode'
            print_string('Preprocessing: pre_proc_norm_acq_parcode()');
            pprocSignals = pre_proc_norm_acq_parcode(sdrParams, caCodeTable);
            
        case  'weak_acq_dbzp'
            print_string('Preprocessing: Calling pre_proc_weak_acq_dbzp()');
            pprocSignals = pre_proc_weak_acq_dbzp(sdrParams, caCodeTable);
            
        otherwise
            error("Acquisition algorithm is not recognized.");
    end
    
    acqAlgoPprocSignals{acqIdx} = pprocSignals;
end

end

