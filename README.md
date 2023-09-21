# RO_prediction_of_PNPA
## Overview
Scripts used to run the retention order prediction presented in the paper:

"Retention order prediction of peptides containing non-progteinogenic amino acids"

_Shohei Nakamukai, Eisuke Hayakawa, Masanori Arita_

## Usage

Retention order prediction of the paper can be reproduced by using the smile_to_prediction.py

1. Ensure you have Python and R installed on your system.

2. Clone this repository or download the source code and data folders.

3. Run the script by executing the following command.

```

import smile_to_prediction
                        
RO_prediction(smile, file)

```

Please enter the SMILES representation of the compound for which you want to predict the holding order in the "smile". In the "file", please enter the absolute path to the folder where you have downloaded the necessary files

For example, when you predict the retention order of α-amanitin and the downloaded files are in "C:/Users/nakamukai/Desktop/", please run the script as follows:

```

import smile_to_prediction

RO_prediction(smile= "CCC(C)[C@@H]1NC(=O)CNC(=O)[C@@H]2Cc3c([nH]c4cc(O)ccc34)S(=O)C[C@H](NC(=O)CNC1=O)C(=O)N[C@@H](CC(=O)N)C(=O)N1C[C@H](O)C[C@H]1C(=O)N[C@@H]([C@@H](C)[C@@H](O)CO)C(=O)N2", file ="C:/Users/nakamukai/Desktop/")

```

The result is displayed as 0 or 1, where 1 indicates a prediction of retention order that comes after the reference compound, and 0 indicates a prediction of retention order that comes before the reference compound.
