
import numpy as np
from rdkit import Chem
#from rdkit.Chem import PandasTools
#
#def pandasDF():
#
#    df = PandasTools.LoadSDF('code/webservice/compute/tsne_cosine_100.sdf', embedProps=True, molColName=None, smilesName='smiles')
#
#    return df

def importSDF():

    suppl = Chem.SDMolSupplier('code/webservice/compute/tsne_cosine_100.sdf')

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
