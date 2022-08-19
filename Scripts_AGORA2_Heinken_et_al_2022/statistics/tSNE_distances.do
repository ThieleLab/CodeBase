*tSNE analyses

*Johannes Hertel

clear
clear mata
clear matrix
set more off
capture log close
set maxvar 40000

cd A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\tSNE_Distances

log using A:\AGORA_2_New\Files_for_Johannes_revision\results\logs\tSNE_distances.log, replace
import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\tSNE_Distances\Bacilli_Reaction_presence_Genus.csv"

loneway yaxisvalue stratifyingfeature
kwallis yaxisvalue, by(stratifyingfeature)
loneway xaxisvalue stratifyingfeature
kwallis xaxisvalue, by(stratifyingfeature)

clear

import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\tSNE_Distances\Gammaproteobacteria_Reaction_presence_Genus.csv"

loneway yaxisvalue stratifyingfeature
kwallis yaxisvalue, by(stratifyingfeature)
loneway xaxisvalue stratifyingfeature
kwallis xaxisvalue, by(stratifyingfeature)

clear

import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\tSNE_Distances\Reaction_presence_Class.csv"

loneway yaxisvalue stratifyingfeature
kwallis yaxisvalue, by(stratifyingfeature)
loneway xaxisvalue stratifyingfeature
kwallis xaxisvalue, by(stratifyingfeature)

clear

import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\tSNE_Distances\Reaction_presence_Family.csv"

loneway yaxisvalue stratifyingfeature
kwallis yaxisvalue, by(stratifyingfeature)
loneway xaxisvalue stratifyingfeature
kwallis xaxisvalue, by(stratifyingfeature)

clear