from flask import Flask, render_template, request, jsonify
from flask_sqlalchemy import SQLAlchemy
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
import base64
import io
import google.generativeai as genai


# Custom property prediction module
from molecule_properties import predict_molecule_properties

app = Flask(__name__)

# Configure Gemini API
GEMINI_API_KEY = ''  # Replace with your actual API key
genai.configure(api_key=GEMINI_API_KEY)

# MySQL Configuration using SQLAlchemy and PyMySQL
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://root:yourpassword@localhost/drug_discovery_db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

db = SQLAlchemy(app)

# Function to get predictive properties from Gemini API
def get_gemini_properties(medicine_name):
    try:
        model = genai.GenerativeModel('gemini-pro')
        prompt = f"""
        Provide a comprehensive analysis of the predictive properties for the medicine: {medicine_name}
        Please include the following detailed information:
        1. Therapeutic Class
        2. Potential Molecular Mechanisms
        3. Predicted Pharmacokinetic Properties
        4. Potential Side Effects
        5. Probable Drug Interactions
        6. Estimated Efficacy Markers
        Format the response as a structured JSON with each category as a key.
        """
        response = model.generate_content(prompt)
        return {
            "therapeutic_class": response.text.split("Therapeutic Class:")[1].split("\n")[0].strip(),
            "molecular_mechanisms": response.text.split("Molecular Mechanisms:")[1].split("\n")[0].strip(),
            "pharmacokinetic_properties": response.text.split("Pharmacokinetic Properties:")[1].split("\n")[0].strip(),
            "potential_side_effects": response.text.split("Side Effects:")[1].split("\n")[0].strip(),
            "drug_interactions": response.text.split("Drug Interactions:")[1].split("\n")[0].strip(),
            "efficacy_markers": response.text.split("Efficacy Markers:")[1].split("\n")[0].strip()
        }
    except Exception as e:
        print(f"Error getting Gemini properties: {e}")
        return None

# Load Chembyl Dataset
def load_dataset():
    try:
        df = pd.read_csv('chembyl_dataset.csv')
        return df
    except Exception as e:
        print(f"Error loading dataset: {e}")
        return None

# Function to generate molecule image
def generate_molecule_image(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol, size=(400, 400))
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        return img_str
    except Exception as e:
        print(f"Error generating molecule image: {e}")
        return None

@app.route('/')
def index():
    dataset = load_dataset()
    if dataset is not None:
        molecules = dataset['SMILES'].tolist()
        medicine_names = dataset['Name'].tolist()
        return render_template('index.html', molecules=molecules, medicine_names=medicine_names)
    return "Error loading dataset"

@app.route('/search_molecule', methods=['POST'])
def search_molecule():
    search_term = request.form.get('search_term')
    dataset = load_dataset()
    result = dataset[
        (dataset['Name'].str.contains(search_term, case=False)) | 
        (dataset['SMILES'] == search_term)
    ]
    if not result.empty:
        molecule = result.iloc[0]
        smiles = molecule['SMILES']
        medicine_name = molecule['Name']
        molecule_image = generate_molecule_image(smiles)
        gemini_properties = get_gemini_properties(medicine_name)
        return render_template('molecule_view.html', 
                               molecule_image=molecule_image,
                               smiles=smiles,
                               medicine_name=medicine_name,
                               medicine_property=molecule['Properties'],
                               gemini_properties=gemini_properties)
    return "Molecule not found", 404

@app.route('/molecule_properties', methods=['POST'])
def molecule_properties():
    smiles = request.form.get('smiles')
    properties = predict_molecule_properties(smiles)
    if properties:
        return jsonify(properties)
    else:
        return jsonify({"error": "Could not predict properties"}), 400

@app.route('/get_gemini_properties', methods=['POST'])
def get_gemini_properties_route():
    medicine_name = request.form.get('medicine_name')
    properties = get_gemini_properties(medicine_name)
    if properties:
        return jsonify(properties)
    else:
        return jsonify({"error": "Could not retrieve properties from Gemini"}), 400

@app.route('/add_compound', methods=['POST'])
def add_compound():
    data = request.json
    compound = data.get('compound')
    try:
        mol = Chem.MolFromSmiles(compound)
        img_str = generate_molecule_image(Chem.MolToSmiles(mol))
        return jsonify({
            'image': img_str,
            'name': Chem.MolToInchiKey(mol)
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 400

if __name__ == '__main__':
    app.run(debug=True)
