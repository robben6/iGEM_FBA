import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import cobra
from cobra.io import read_sbml_model
ROOT_DIR = Path(r"C:\Users\DELL\OneDrive - UBC\UBC\Design Team\iGEM")
RAW_DATA_DIR = ROOT_DIR / "Project"/ "FBA modeling"
model = read_sbml_model(RAW_DATA_DIR / "UTEX2973_S2.xml")

# store the knockout results
result = []

solution = model.optimize()
result.append({"Condition": "Complete model", "Objective": solution.objective_value})
print(solution)

##Single gene knockout
print('complete model: ', model.optimize())

with model:
    model.genes.SYNPCC7942_0201.knock_out()
    print('O2_gene1 knocked out: ', model.optimize())

    model.genes.SYNPCC7942_2561.knock_out()
    print('O2_gene2 knocked out: ', model.optimize())

    model.genes.SYNPCC7942_0978.knock_out()  # oxidized
    print('O2_gene3 knocked out: ', model.optimize())

    model.genes.SYNPCC7942_0713.knock_out()  # oxidized
    print('O2_gene4 knocked out: ', model.optimize())

##Single reaction knockout
print('complete model: ', model.optimize())

with model:  # knockout the oxygen exchange reaction
    # model.reactions.CYO1b_syn.knock_out()
    # print('CYO1b_syn knocked out: ', model.optimize())

    # model.reactions.Desat9.knock_out()
    # print('Desat9 knocked out: ', model.optimize())

    model.reactions.EX_O2.knock_out()
    solution = model.optimize()
    result.append({"Condition": "O2 exchange KO", "Objective": solution.objective_value})
    print('oxygen exchange knocked out: ', solution)

with model:  # knockout the oxygen transport

    model.reactions.rxn09031.knock_out()
    solution = model.optimize()
    result.append({"Condition": "O2 transport KO", "Objective": solution.objective_value})
    print('oxygen transport knocked out: ', solution)

## Visualization
# Convert to DataFrame
df = pd.DataFrame(result)
# Bar plot
plt.bar(df["Condition"], df["Objective"])
plt.ylabel("Biomass Objective Value")
plt.title("Effect of Knockouts on Growth")
plt.xticks(rotation=20)
plt.show()

##FVA 
BG_11_MEDIA_EXRXNS_IDS = {
    "EX_CO2",
    "EX_PHO1",
    "EX_PHO2",
    "EX_Citrate",
    "EX_NH3",
    "EX_Mn2_",
    "EX_Zinc",
    "EX_Cu2_",
    "EX_Molybdate",
    "EX_Co2_",
    "EX_Nitrate",
    "EX_Phosphate",
    "EX_Sulfate",
    "EX_Fe3_",
    "EX_H2CO3",
    "EX_Calcium"
}
for rxn_id in BG_11_MEDIA_EXRXNS_IDS:
    if rxn_id not in model.reactions:
        raise ValueError(f"Reaction {rxn_id} not found in the model.")
    else:
        print(model.reactions.get_by_id(rxn_id))

