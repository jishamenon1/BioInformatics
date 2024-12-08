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
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.preprocessing import StandardScaler
import random

# Load the trained model and scaler
model = joblib.load('if_predictor_model.pkl')
scaler = joblib.load('scaler.pkl')

# Function to convert SMILES to descriptors
def smiles_to_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return [Descriptors.MolWt(mol), Descriptors.MolLogP(mol)]
    else:
        return [np.nan, np.nan]

# Add custom styling with background animation and circular balloon effect
st.markdown(
    """
    <style>
    h1 {
        color: #00E676;
        font-family: 'Arial', sans-serif;
        text-align: center;
    }

    h2 {
        color: #00E676;
        font-family: 'Arial', sans-serif;
        text-align: center;
    }

    .section-title {
        font-size: 18px;
        color: #9457eb;
        font-weight: bold;
        margin-top: 10px;
    }

    .prediction {
        color: #00E676;
        text-align: center;
        font-size: 24px;
        margin-top: 0;  /* Remove top margin */
        padding: 10px 0;  /* Add small padding to center the text better */
    }

    .error {
        color: #FF5252;
        font-weight: bold;
    }

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
    """,
    unsafe_allow_html=True
)

st.markdown('<div class="main-container">', unsafe_allow_html=True)

# App title and description
st.markdown('<h2>Predict IF Value from SMILES and Concentrations Using ML</h2>', unsafe_allow_html=True)
st.write("🔬 Provide the SMILES  and their concentrations to predict the Imprinting Factor (IF).")

# Input section
st.markdown('<div class="section-title">Input Details</div>', unsafe_allow_html=True)

template_smiles = st.text_input("🧬 Template SMILES")
template_conc = st.number_input("💧 Template Concentration (mmol)", min_value=0.0, step=0.01)

crosslinker_smiles = st.text_input("🔗 Crosslinker SMILES")
crosslinker_conc = st.number_input("💧 Crosslinker Concentration (mmol)", min_value=0.0, step=0.01)

monomer_smiles = st.text_input("🧪 Monomer SMILES")
monomer_conc = st.number_input("💧 Monomer Concentration (mmol)", min_value=0.0, step=0.01)

# Prediction button with balloon animation
if st.button("📊 Predict IF Value"):
    # Only trigger balloons if all SMILES are provided
    if template_smiles and crosslinker_smiles and monomer_smiles:
        for i in range(5):  # Show 5 colorful balloons for a more exciting effect
            st.markdown(f'<div class="balloon balloon-{i+1} show">🎈</div>', unsafe_allow_html=True)

        # Convert SMILES to descriptors
        template_desc = smiles_to_descriptors(template_smiles)
        crosslinker_desc = smiles_to_descriptors(crosslinker_smiles)
        monomer_desc = smiles_to_descriptors(monomer_smiles)

        if np.isnan(template_desc).any() or np.isnan(crosslinker_desc).any() or np.isnan(monomer_desc).any():
            st.error("❌ Invalid SMILES detected. Please check your inputs.")
        else:
            # Prepare input data (only 8 features)
            input_data = np.array([
                template_desc[0], template_desc[1], template_conc,
                crosslinker_desc[0], crosslinker_desc[1], crosslinker_conc,
                monomer_desc[0], monomer_desc[1]
            ]).reshape(1, -1)  # Ensure 8 features only

            # Scale and predict
            try:
                input_scaled = scaler.transform(input_data)
                predicted_if = model.predict(input_scaled)[0]

                # Display the prediction without a border or background
                st.markdown(f"<div class='prediction'>🌟 Predicted IF Value: {predicted_if:.2f}</div>", unsafe_allow_html=True)
                
                # Button for new prediction
                if st.button("🔄 New Prediction"):
                    # Clear all inputs and reset the page for a new prediction
                    st.experimental_rerun()

            except ValueError as e:
                st.error(f"❌ Error during prediction: {e}")
    else:
        st.error("❌ Please provide valid SMILES strings for Template, Crosslinker, and Monomer.")

st.markdown('</div>', unsafe_allow_html=True)