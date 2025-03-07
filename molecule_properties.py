import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, AllChem

def predict_molecule_properties(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        
        # Comprehensive molecular property calculation
        properties = {
            # Basic Molecular Descriptors
            "Molecular Weight": round(Descriptors.ExactMolWt(mol), 2),
            "Molecular Formula": AllChem.CalcMolFormula(mol),
            "LogP": round(Crippen.MolLogP(mol), 2),
            
            # Structural Properties
            "Hydrogen Bond Donors": Descriptors.NumHDonors(mol),
            "Hydrogen Bond Acceptors": Descriptors.NumHAcceptors(mol),
            "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
            "Topological Polar Surface Area": round(Descriptors.TPSA(mol), 2),
            
            # Toxicity Prediction Indicators
            "Lipinski's Rule of Five": check_lipinski_rule(mol),
            "Potential Toxicity Risk": assess_toxicity_risk(mol),
            
            # Pharmacological Indicators
            "Drug-likeness Score": calculate_drug_likeness(mol)
        }
        
        return properties
    
    except Exception as e:
        print(f"Error in property prediction: {e}")
        return None

def check_lipinski_rule(mol):
    """
    Lipinski's Rule of Five: Criteria for oral bioavailability
    - Molecular weight < 500 Da
    - LogP < 5
    - Hydrogen bond donors ≤ 5
    - Hydrogen bond acceptors ≤ 10
    """
    mol_weight = Descriptors.ExactMolWt(mol)
    log_p = Crippen.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    
    violations = 0
    if mol_weight > 500: violations += 1
    if log_p > 5: violations += 1
    if h_donors > 5: violations += 1
    if h_acceptors > 10: violations += 1
    
    return f"{4 - violations}/4 Rules Passed"

def assess_toxicity_risk(mol):
    """
    Basic toxicity risk assessment based on molecular descriptors
    """
    mol_weight = Descriptors.ExactMolWt(mol)
    log_p = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    
    # Very basic risk scoring (this is a simplified mock)
    risk_score = 0
    if mol_weight > 500: risk_score += 1
    if log_p > 5: risk_score += 1
    if tpsa > 140: risk_score += 1
    
    risk_levels = {
        0: "Low Risk",
        1: "Moderate Risk",
        2: "High Risk",
        3: "Very High Risk"
    }
    
    return risk_levels.get(risk_score, "Unknown Risk")

def calculate_drug_likeness(mol):
    """
    Calculate a simple drug-likeness score
    """
    # Combination of various descriptors
    mol_weight = Descriptors.ExactMolWt(mol)
    log_p = Crippen.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    
    # Simplified scoring (this is a mock implementation)
    score = 0
    if 200 < mol_weight < 500: score += 1
    if -0.4 < log_p < 5.6: score += 1
    if h_donors <= 5: score += 1
    if h_acceptors <= 10: score += 1
    if rotatable_bonds <= 10: score += 1
    
    return f"{score}/5 Drug-likeness"