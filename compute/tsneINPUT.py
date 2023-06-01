import numpy as np
import pandas as pd
from rdkit import Chem

def importCSV(tsneDirName, tsneFileName):

    from rdkit.Chem import Descriptors
    from rdkit.Chem import rdMolDescriptors
    #from rdkit.Chem.Draw import rdMolDraw2D
    #from rdkit.Chem import Draw
    #from rdkit.Chem.Draw import MolToFile

    df = pd.read_csv("code/webservice/compute/data/tsne/"+tsneDirName+"/"+tsneFileName+".csv")

    graphX = []
    graphY = []
    graphZ = []
    smilesList = []
    molForList = []
    molWeightList = []

    allvariables = []
    for var in df:
        allvariables.append(str(var))

    for mol in range(len(df)):
        graphX.append(float(df.at[mol, allvariables[1]]))
        graphY.append(float(df.at[mol, allvariables[2]]))
        graphZ.append(float(df.at[mol, allvariables[4]]))
        smilesList.append(str(df.at[mol, allvariables[3]]))

        rdkitMol =  Chem.MolFromSmiles(str(df.at[mol, allvariables[3]]))

        molForList.append(str(rdMolDescriptors.CalcMolFormula(rdkitMol)))
        molWeightList.append(str(Descriptors.MolWt(rdkitMol)))

        #img = Chem.Draw.MolToImage(rdkitMol)


    #    #print(df.at[mol, 'Unnamed: 0'],df.at[mol, 'x'], df.at[mol, 'y'], df.at[mol, 'Smiles'], df.at[mol, 'MAE'])

    molSmiles = str(smilesList).replace("'","")
    molFormula = str(molForList).replace("'","")
    molWeight = str(molWeightList).replace("'","")

    return graphX, graphY, graphZ, molFormula, molSmiles, molWeight

def importSDF():

    suppl = Chem.SDMolSupplier('code/webservice/compute/data/tsne_cosine_100.sdf')

    graphX = []
    graphY = []
    graphZ = []
    molForList = []
    smilesList = []
    molWeightList = []
    
    for mol in suppl:
        graphX.append(float(mol.GetProp('x')))
        graphY.append(float(mol.GetProp('y')))
        graphZ.append(float(mol.GetProp('Label')))
        molForList.append(str(mol.GetProp('Molecular Formula')))
        smilesList.append(str(mol.GetProp('Smiles')))
        molWeightList.append(str(mol.GetProp('Total Molweight')))



    #remove ' from string (not know how to properly do it in js)
    molFormula = str(molForList).replace("'","")
    molSmiles = str(smilesList).replace("'","")
    molWeight = str(molWeightList).replace("'","")

    return graphX, graphY, graphZ, molFormula, molSmiles, molWeight

def gaussPoints(minX, maxX, npoints, a, b, c):

    gaussx=[]
    gaussy=[]

    delta = (maxX - minX) / npoints

    for x in np.arange(minX, maxX+delta, delta):

        y = a * np.exp( -np.power(x-b,2) / (2.0*np.power(c,2))  )

        gaussx.append(x)
        gaussy.append(y)

    return gaussx, gaussy
