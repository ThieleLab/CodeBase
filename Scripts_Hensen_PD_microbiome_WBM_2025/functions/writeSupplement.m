function writeSupplement(suplTable, description, folder, suplTable2)
% Helper function for populating the supplementary materials excel file

if nargin<4
    suplTable2 = '';
end

% Define file path
filePath = fullfile(folder,'Supplementary_tables.xlsx');

% Find the current sheets if the supplementary table file exists
SM = isfile(filePath);
if SM == true
    sheets = sheetnames(filePath); % Find the current sheet names
    maxSheetNum = max(str2double(erase(sheets,'Table_S'))); % Find the largest sheet number
    if ~isnan(maxSheetNum) % Double check
        sheetNum = char(string(maxSheetNum+1)); % Define the current sheet number as the max + 1
    else
        sheetNum = '1'; % Set sheet number to 1
    end
else
    sheetNum = '1'; % Set sheet number to 1
end

% Create table header:
sheetName = append('Table_S',sheetNum); % Process sheet name
tableHeader = append(sheetName,": ",string(description{1}));

% Add details to table header
if ~isempty(description{2}) % Only add extra line if a description is given
    details = append("Description: ",string(description{2}));
else
    details = string(description{2});
end

tableHeader = [tableHeader; details];


% Create of update index sheet
if SM == false
    tableHeader = ["INDEX"; tableHeader]; % Add Table header if not present already
end

% Write index table to file, but remove the table details
writematrix(tableHeader( 1:(end-1) ), filePath,'Sheet','Index','WriteMode','append') 

% Remove Index table header
tableHeader(matches(tableHeader,"INDEX"))=[];

% Create new sheet for current supplementary table and write table
% description
writematrix(tableHeader, filePath,'Sheet',sheetName,'WriteMode','overwritesheet')

% Append supplementary table to excel sheet
writetable(suplTable,filePath,'Sheet',sheetName,'Range','A4')

if ~isempty(suplTable2)
    % Place the second table after the second empty column in excel sheet after the main table:
    colNum = width(suplTable)+3;
    colLetter = char(colNum + 64); % Convert to letter in alphabet using ASCII codes
    rangeStart = [colLetter,'4']; % Find excel cell to start
    
    % Append the second supplementary table to excel sheet
    writetable(suplTable2,filePath,'Sheet',sheetName,'Range',rangeStart)
end

disp(append('Saved ',sheetName, ' to Supplementary_tables.xlsx'))
end