
import streamlit as st
from stmol import showmol
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp


# Custom CSS for styling
st.markdown("""
    <style>
        body {
            font-family: Arial, sans-serif;
            background-color: #f9f9f9;
        }
        .main-title {
            text-align: center;
            color: #4CAF50;
            font-size: 3rem;
            font-weight: bold;
            margin-bottom: 1rem;
        }
        .sub-title {
            text-align: center;
            color: #555;
            font-size: 1.2rem;
            margin-bottom: 2rem;
        }
       
        .section h2 {
            color: #4CAF50;
            font-size: 1.8rem;
        }
        .section hr {
            border: none;
            height: 1px;
            background-color: #ddd;
            margin: 1rem 0;
        }
        .footer {
            text-align: center;
            color: #aaa;
            font-size: 0.9rem;
            margin-top: 2rem;
        }
    </style>
""", unsafe_allow_html=True)

# Title Section
st.markdown('<div class="main-title">3D Molecule Viewer</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Visualize molecular structures from SMILES </div>', unsafe_allow_html=True)

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


# Helper Functions
def makeblock(smi):
    """Convert SMILES to a 3D MolBlock."""
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def render_mol(xyz):
    """Render a molecule in 3D using py3Dmol."""
    xyzview = py3Dmol.view()
    xyzview.addModel(xyz, 'mol')
    xyzview.setStyle({'stick': {}})
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    # Adjust the size of the visualizer
    showmol(xyzview, height=600, width=700)

# Section 1: SMILES Input and 3D Model Rendering
with st.container():
    st.markdown('<div class="section">', unsafe_allow_html=True)
    st.subheader("🔬 Generate 3D Model from SMILES")

    # Input field with placeholder that disappears on typing
    compound_smiles = st.text_input(
        "Enter SMILES", 
        placeholder="Example: CC(C)(C)C1=CC=C(C=C1)O (click to type)",  # Placeholder text with example
        help="Input a SMILES string to visualize its 3D structure."
    )

  # Button to trigger rendering
    if st.button("Generate 3D Model"):
        if compound_smiles:  # Check if SMILES input is provided
            blk = makeblock(compound_smiles)
            render_mol(blk)
        else:
            st.warning("Please enter a valid SMILES string before clicking the button.")
# Footer


