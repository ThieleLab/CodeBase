function metadata = processMetadata(metadataPath)
% processMetadata checks if the required variables, ID and sex, are
% included in the metadata. If the ID column is called, sample, name,
% sample_id, or samples_name, the variable is renamed to ID. If the
% variable sex encodes sex information in the form, m/f, the entries are
% translated to "male" and "female". An error is thrown if the metadata
% contains names or formats outside the allowed formats.
%
% INPUT
% metadataPath          Path to the metadata file
% 
% OUTPUT
% updatedMetadataPath   Path to the processed metadata file
%
% AUTHOR: Tim Hensen, July 2024

% Read the metadata file
metadata = readtable(metadataPath,'VariableNamingRule','preserve');

% Make variable names work with matlab
metadata.Properties.VariableNames = matlab.lang.makeValidName(metadata.Properties.VariableNames);

% Check if the sample IDs and Sex information are included
%metadata.Properties.VariableNames = lower(metadata.Properties.VariableNames);
varNames = metadata.Properties.VariableNames;
IDcol = lower(metadata.Properties.VariableNames{1});
% Check if the sample names are correct and correct the variable name if
% necessary.
switch IDcol
    case 'id'
        metadata = renamevars(metadata,'id','ID');
    case 'sample'
        metadata = renamevars(metadata,'sample','ID');
    case 'name'
        metadata = renamevars(metadata,'name','ID');
    case 'sample_id'
        metadata = renamevars(metadata,'sample_id','ID');
    case 'sample_name'
        metadata = renamevars(metadata,'sample_name','ID');
    otherwise
        error('Cannot find sample IDs. Make sure that the sample IDs are in the first column and are named ID')
end

% Check if a variable exists with sex information
if any(matches(lower(varNames),'gender'))
    metadata = renamevars(metadata,'gender','Sex');
elseif ~matches(lower(varNames),'sex')
    error('Cannot find sample Sex variable. Make sure that the sample Sex is in the metadata')
end

% Ensure the IDs and Sex info are string arrays
metadata.ID = string(metadata.ID);
metadata.Sex = string(metadata.Sex);

% Check if male and/or female samples can be found
metadata.Sex = lower(metadata.Sex);
if any(matches(metadata.Sex,{'f','m'}))
    % Convert F and M to female and female
    metadata.Sex(matches(metadata.Sex,'f')) = "female";
    metadata.Sex(matches(metadata.Sex,'m')) = "male";
end

if ~any(matches(metadata.Sex,{'female','male'}))
    error('Please ensure that sex information is encoded as "female" and "male".')
end

end
