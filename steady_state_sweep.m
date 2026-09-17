%% Run steady-state search over all parameter combinations

clear
clc

% Parameter values
%pVals      = 0:0.1:1;
%gammaVals  = 10.^(-6:0.5:-1);
%RetaVals   = [0 10 100 1000];
%filterVals = [1 2 3];
 
%paramter vals for test
pVals      = 0:.5:1;
gammaVals  = 10.^(-6:3:-1);
RetaVals   = [0 10];
filterVals = [1 2];

% Create list of all parameter combinations
nP      = length(pVals);
nGamma  = length(gammaVals);
nReta   = length(RetaVals);
nFilter = length(filterVals);

nComb = nP * nGamma * nReta * nFilter;

% Preallocate arrays for parameter combinations
pList     = zeros(nComb,1);
gammaList = zeros(nComb,1);
RetaList  = zeros(nComb,1);
filterList = zeros(nComb,1);

counter = 0;

for ip = 1:nP
    for ig = 1:nGamma
        for ir = 1:nReta
            for ifilt = 1:nFilter

                counter = counter + 1;

                pList(counter)      = pVals(ip);
                gammaList(counter)  = gammaVals(ig);
                RetaList(counter)   = RetaVals(ir);
                filterList(counter) = filterVals(ifilt);

            end
        end
    end
end

% Run all combinations

resultsCell = cell(nComb,1);

parfor i = 1:nComb

    resultsCell{i} = finding_steady_state( ...
        pList(i), ...
        gammaList(i), ...
        RetaList(i), ...
        filterList(i));

end

% Combine all results into one table
results = vertcat(resultsCell{:});

% Save final table
writetable(results, 'steady_state_results.csv');

save('steady_state_results.mat', 'results');
