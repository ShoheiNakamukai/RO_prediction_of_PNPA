# RO_prediction_of_PNPA
## Overview
Scripts used to run the retention order prediction presented in the paper:

"Retention order prediction of peptides containing non-progteinogenic amino acids"

_Shohei Nakamukai, Eisuke Hayakawa, Masanori Arita_

## Usage

Retention order prediction of the paper can be reproduced by using the smile_to_prediction.py

1. Ensure you have Python and R installed on your system.

2. Clone this repository or download the source code and data folders.

3. Prepare the all scripts from [here](https://version.aalto.fi/gitlab/bache1/retention_order_prediction/-/tree/master/src?ref_type=heads). Please place all of them in the folder "RO_prediction_of_PNPA/src". 

4. Run the script by executing the following command.

```

import smile_to_prediction
                        
smile_to_prediction.RO_prediction(smile, file_name)

```

Please input the "smile" with the SMILES representation of the molecule whose retention order you want to predict. In the "file", please enter the absolute path to the folder where you have downloaded the necessary files

The result is displayed as 0 or 1, where 1 indicates a prediction of retention order that comes after the reference compound, and 0 indicates a prediction of retention order that comes before the reference compound.

## Example

For example, when you predict the retention order of α-amanitin and the downloaded files are in "C:/Users/nakamukai/Desktop/RO_prediction_of_PNPA/", please run the script as follows:

```

import smile_to_prediction

smile_to_prediction.RO_prediction(smile= "CCC(C)[C@@H]1NC(=O)CNC(=O)[C@@H]2Cc3c([nH]c4cc(O)ccc34)S(=O)C[C@H](NC(=O)CNC1=O)C(=O)N[C@@H](CC(=O)N)C(=O)N1C[C@H](O)C[C@H]1C(=O)N[C@@H]([C@@H](C)[C@@H](O)CO)C(=O)N2", file_name ="C:/Users/nakamukai/Desktop/RO_prediction_of_PNPA/")

```

You can obtain the result as follows:

```

2-Propenyl glucosinolate: 1
p-Hydroxybenzyl glucosinolate: 1
Cyanidin-3-O-sambubioside-5-O-glucoside: 1
Carbazochrome sulfonate: 1
Sinomenine: 1
(-)-Epicatechin: 0
7,8-Dihydroxycoumarin: 0
Isovitexin(4): 0
3-Hydroxycinnamic acid: 0
Coniferyl aldehyde: 0
3,4-Dimethoxycinnamic acid: 0
Quercetin: 0
4-Methoxycinnamic acid: 0
2-Methoxycinnamic acid: 0
isoliquiritigenin: 0
6-Hydroxyflavanone: 0
2'-Hydroxyflavanone: 0
Atractylenolide III: 0
trans-pterostilbene: 0
Triacetyl resveratrol: 0
Magnolol: 0
Corosolic acid: 0

```

The result means that α-amanitin will elute between 3-Hydroxycinnamic acid and Coniferyl aldehyde in the reversed phase liquid chromatography.

## Reference

Scripts in the src folder are cited from []
The SMARTS pattern in MACCS_list_166.csv can be downloaded from [here](https://github.com/cdk/cdk/blob/4004eb64fd7e94a0da674ae2c0eedba79fda425f/descriptor/fingerprint/src/main/resources/org/openscience/cdk/fingerprint/data/SMARTS_countable_MACCS_keys.txt)
