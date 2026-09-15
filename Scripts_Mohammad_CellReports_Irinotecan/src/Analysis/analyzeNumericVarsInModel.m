function [numericInfoCellArray] = analyzeNumericVarsInModel(model)
    % Check if the model structure is empty
    if isempty(model)
        disp('Model structure is empty.');
        numericInfoCellArray = [];
        return;
    end
    
    % Display a message indicating the presence of the model structure
    disp('Model structure is not empty.');

    % Initialize a cell array to store numeric variables and their max digit lengths
    numericInfoCellArray = {};

    % Iterate through all fields in the model structure
    fieldNames = fieldnames(model);
    for i = 1:numel(fieldNames)
        varName = fieldNames{i};
        varValue = model.(varName);
        numericInfo = analyzeVariable(varName, varValue);

        % Store numeric variables and their max digit lengths in the cell array
        if ~isempty(numericInfo)
            numericInfoCellArray = [numericInfoCellArray; numericInfo];
        end
    end

    % Display information about numeric variables and their max digit lengths
    displayNumericInfo(numericInfoCellArray);
end

function numericInfo = analyzeVariable(varName, varValue)
    % Check if the variable is present
    if isempty(varValue)
        fprintf('%s: Variable is empty.\n', varName);
        numericInfo = [];
        return;
    end

    % Determine if the variable is sparse
    isSparse = issparse(varValue);

    % Calculate the maximum digit lengths for a variable
    if isnumeric(varValue)
        if isSparse
            maxDigitLength = max(arrayfun(@numDigits, nonzeros(varValue(:))));
        else
            maxDigitLength = max(arrayfun(@numDigits, varValue(:)));
        end
        fprintf('%s: Maximum digit length is %d.\n', varName, maxDigitLength);

        % Store numeric variable name and its max digit length in a cell array
        numericInfo = {varName, maxDigitLength};
    else
        fprintf('%s: Not a numeric variable.\n', varName);
        numericInfo = [];
    end
end

function digits = numDigits(value)
    % Helper function to calculate the number of digits in a numeric value
    digits = numel(num2str(value));
end

function displayNumericInfo(numericInfoCellArray)
    % Display information about numeric variables and their max digit lengths
    fprintf('\nNumeric Variables Information:\n');
    if isempty(numericInfoCellArray)
        disp('No numeric variables found.');
    else
        for i = 1:size(numericInfoCellArray, 1)
            varName = numericInfoCellArray{i, 1};
            maxDigitLength = numericInfoCellArray{i, 2};
            fprintf('%s: Maximum digit length is %d.\n', varName, maxDigitLength);
        end
    end
end
