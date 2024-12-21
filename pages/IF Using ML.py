import streamlit as st
import pubchempy as pcp

# Function to fetch SMILES from PubChem using a compound name
def fetch_smiles_from_pubchem(compound_name):
    try:
        # Search for the compound by name
        compounds = pcp.get_compounds(compound_name, 'name', timeout=10)
        
        # Check how many compounds were found
        num_compounds = len(compounds)
        
        if num_compounds > 0:
            compound = compounds[0]  # Get the first result
            smiles = compound.isomeric_smiles  # Get the isomeric SMILES
            return smiles, num_compounds, compound
        else:
            return None, num_compounds, None
    except Exception as e:
        return f"Error fetching data for {compound_name}: {str(e)}", 0, None

# Streamlit UI setup
st.markdown(
    """
    <style>
    /* Main Background Styling */
    .main {{
        background-color: #4CAF50;  /* Green background color */
        color: #FFFFFF;  /* White text color */
        padding: 2rem;
        border-radius: 10px;
        text-align: center;
        font-family: 'Arial', sans-serif;
    }}
    
    h1, h2 {{
        color: #FFFFFF;  /* White color for headings */
        font-family: 'Arial', sans-serif;
        font-weight: bold;
    }}
    
    .stButton button {{
        background-color: #FFC107;  /* Yellow background for button */
        color: white;
        border-radius: 5px;
        font-size: 18px;
        width: 200px;
    }}
    .stButton button:hover {{
        background-color: #FF9800;  /* Orange on hover */
    }}

    .stTextInput input {{
        font-size: 18px;
        border-radius: 5px;
        padding: 0.5rem;
        width: 100%;
        border: 1px solid #ccc;
    }}

    .stText {{
        font-size: 18px;
    }}
    </style>
    """, unsafe_allow_html=True)

# Main content div
st.markdown('<div class="main">', unsafe_allow_html=True)
# Input from user for compound name
with st.form(key='compound_form'):
    compound_name = st.text_input("🔍 Enter compound name:")
    submit_button = st.form_submit_button("📥 Get SMILES")

# Button to fetch the SMILES
if submit_button:
    if compound_name:
        with st.spinner('Fetching SMILES ...'):
            smiles, num_compounds, compound = fetch_smiles_from_pubchem(compound_name)
            
            if smiles:
                st.write(f"**SMILES for {compound_name}:**")
                st.text(smiles)  # Display the SMILES string
               
            else:
                st.write(f"No compounds found for {compound_name}.")
    else:
        st.write("Please enter a valid compound name.")

import streamlit as st
import pubchempy as pcp
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import py3Dmol
from stmol import showmol
from sklearn.preprocessing import StandardScaler

# Initialize session state variables
for key in [
    "template_smiles",
    "template_conc",
    "crosslinker_smiles",
    "crosslinker_conc",
    "monomer_smiles",
    "monomer_conc",
    "prediction_complete",
]:
    if key not in st.session_state:
        st.session_state[key] = "" if "smiles" in key else 0.0 if "conc" in key else False

# CSS Styling
st.markdown("""
    <style>
    h1, h2 { text-align: center; color: #00E676; }
    .prediction { color: #00E676; text-align: center; font-size: 24px; display: flex; align-items: center; justify-content: center; }
    .error { color: #FF5252; font-weight: bold; }
    .structure-row { display: flex; justify-content: center; gap: 30px; margin-top: 20px; }

    /* Balloon Animation */
    @keyframes balloonFloat {
        0% {
            transform: scale(0) translate(0, 0);
            opacity: 1;
        }
        25% {
            transform: scale(1) translate(80px, -120px);
            opacity: 0.9;
        }
        50% {
            transform: scale(1.2) translate(120px, -200px);
            opacity: 0.8;
        }
        75% {
            transform: scale(1) translate(-100px, -250px);
            opacity: 0.6;
        }
        100% {
            transform: scale(1.5) translate(-150px, -300px);
            opacity: 0;
        }
    }

    .balloon {
        font-size: 30px;
        position: absolute;
        animation: balloonFloat 3s ease-out forwards;
        display: none;
    }

    .balloon.show {
        display: block;
    }

    /* Balloon color classes */
    .balloon-1 { top: 50%; left: 50%; animation-delay: 0.1s; transform-origin: -80px -80px; color: #FF6347; } /* Red */
    .balloon-2 { top: 50%; left: 50%; animation-delay: 0.2s; transform-origin: 100px 100px; color: #FFD700; } /* Yellow */
    .balloon-3 { top: 50%; left: 50%; animation-delay: 0.3s; transform-origin: -120px 120px; color: #32CD32; } /* Green */
    .balloon-4 { top: 50%; left: 50%; animation-delay: 0.4s; transform-origin: 140px -140px; color: #1E90FF; } /* Blue */
    .balloon-5 { top: 50%; left: 50%; animation-delay: 0.5s; transform-origin: 160px -160px; color: #FF1493; } /* Pink */

    
    </style>
""", unsafe_allow_html=True)

# Functions
def smiles_to_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return [Descriptors.MolWt(mol), Descriptors.MolLogP(mol)]
    return [np.nan, np.nan]

def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    return Chem.MolToMolBlock(mol)

def render_mol(mol_block):
    view = py3Dmol.view()
    view.addModel(mol_block, "mol")
    view.setStyle({"stick": {}})
    view.setBackgroundColor("white")
    view.zoomTo()
    return view

# Title
st.markdown('<h2>Predict IF Value and Visualize Molecular Structures</h2>', unsafe_allow_html=True)

# Inputs
st.session_state["template_smiles"] = st.text_input("🧬 Template SMILES", st.session_state["template_smiles"])
st.session_state["template_conc"] = st.number_input("💧 Template Concentration (mmol)", min_value=0.0, step=0.01, value=st.session_state["template_conc"])
st.session_state["crosslinker_smiles"] = st.text_input("🔗 Crosslinker SMILES", st.session_state["crosslinker_smiles"])
st.session_state["crosslinker_conc"] = st.number_input("💧 Crosslinker Concentration (mmol)", min_value=0.0, step=0.01, value=st.session_state["crosslinker_conc"])
st.session_state["monomer_smiles"] = st.text_input("🧪 Monomer SMILES", st.session_state["monomer_smiles"])
st.session_state["monomer_conc"] = st.number_input("💧 Monomer Concentration (mmol)", min_value=0.0, step=0.01, value=st.session_state["monomer_conc"])

# Prediction and 3D Visualization
if st.button("📊 Predict IF Value"):
    # Show balloons when the button is clicked
    st.session_state["show_balloons"] = True
    
    if st.session_state["template_smiles"] and st.session_state["crosslinker_smiles"] and st.session_state["monomer_smiles"]:
        template_desc = smiles_to_descriptors(st.session_state["template_smiles"])
        crosslinker_desc = smiles_to_descriptors(st.session_state["crosslinker_smiles"])
        monomer_desc = smiles_to_descriptors(st.session_state["monomer_smiles"])

        if np.isnan(template_desc).any() or np.isnan(crosslinker_desc).any() or np.isnan(monomer_desc).any():
            st.error("❌ Invalid SMILES detected. Please check your inputs.")
        else:
            input_data = np.array([
                template_desc[0], template_desc[1], st.session_state["template_conc"],
                crosslinker_desc[0], crosslinker_desc[1], st.session_state["crosslinker_conc"],
                monomer_desc[0], monomer_desc[1]
            ]).reshape(1, -1)

            try:
                model = joblib.load("if_predictor_model.pkl")
                scaler = joblib.load("scaler.pkl")
                input_scaled = scaler.transform(input_data)
                predicted_if = model.predict(input_scaled)[0]

                # Show prediction result and balloon animation
                with st.container():
                    # Show balloon animation next to the predicted value
                    st.markdown(f"<div class='prediction'>🌟 Predicted IF Value: {predicted_if:.2f} <div class='balloon balloon-1 show'>🎈</div><div class='balloon balloon-2 show'>🎈</div><div class='balloon balloon-3 show'>🎈</div><div class='balloon balloon-4 show'>🎈</div><div class='balloon balloon-5 show'>🎈</div></div>", unsafe_allow_html=True)
                    
                    # Show 3D structure in a white box with increased width
                    st.markdown('<div class="structure-box">', unsafe_allow_html=True)
                    for name, smiles in [("Template", st.session_state["template_smiles"]), 
                                         ("Crosslinker", st.session_state["crosslinker_smiles"]), 
                                         ("Monomer", st.session_state["monomer_smiles"])]:
                        mol_block = makeblock(smiles)
                        st.write(f"### {name}")
                        showmol(render_mol(mol_block), height=400, width=700)  # Increased width here
                    st.markdown('</div>', unsafe_allow_html=True)

            except Exception as e:
                st.error(f"❌ Error during prediction: {e}")
    else:
        st.error("❌ Please provide valid SMILES strings for Template, Crosslinker, and Monomer.")


