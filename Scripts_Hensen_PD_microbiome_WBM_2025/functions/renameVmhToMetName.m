function metaboliteList = renameVmhToMetName(vmhList, reverse)

if nargin <2
    reverse = false;
end

% Rename VMH reactions to metabolites
vmhList = string(vmhList);
vmhList(matches(vmhList,"DM_leuleu[bc]")) = "Leucylleucine";
vmhList(matches(vmhList,"DM_leu_L[bc]")) = "L-leucine";
vmhList(matches(vmhList,"DM_but[bc]")) = "Butyrate";
vmhList(matches(vmhList,"DM_pnto_R[bc]")) = "Pantothenate";
vmhList(matches(vmhList,"DM_nac[bc]")) = "Nicotinic acid";
vmhList(matches(vmhList,"DM_ttdca[bc]")) = "Myristic acid";
metaboliteList = cellstr(vmhList);

if reverse == true
    % Rename metabolites to VMH reactions
    vmhList = string(vmhList);
    vmhList(matches(vmhList,"Leucylleucine")) = "DM_leuleu[bc]";
    vmhList(matches(vmhList,"L-leucine")) = "DM_leu_L[bc]";
    vmhList(matches(vmhList,"Butyrate")) = "DM_but[bc]";
    vmhList(matches(vmhList,"Pantothenate")) = "DM_pnto_R[bc]";
    vmhList(matches(vmhList,"Nicotinic acid")) = "DM_nac[bc]";
    vmhList(matches(vmhList,"Myristic acid")) = "DM_ttdca[bc]";
    metaboliteList = cellstr(vmhList);

end