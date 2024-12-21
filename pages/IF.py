import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# Load dataset
dataset = pd.read_csv('your_dataset.csv')  # Replace with your dataset file path
dataset.columns = dataset.columns.str.strip()  # Clean column names

# Custom CSS for styling
st.markdown("""
    <style>
        body {
            background-color: #f9f9fb;
            color: #333;
            font-family: 'Arial', sans-serif;
        }
        .header {
            font-size: 2.5rem;
            font-weight: bold;
            text-align: center;
            color: #4CAF50;
            margin-bottom: 1rem;
        }
        .sub-header {
            text-align: center;
            color: #555;
            font-size: 1.2rem;
            margin-bottom: 2rem;
        }
        .if-value {
            font-size: 1.5rem;
            font-weight: bold;
            color: #FF6347;  
            background-color: #FFE4E1; 
            padding: 7px;
            border-radius: 4px;
            text-align: center;
            margin-bottom: 0.5rem;
        }
        .details-title {
            font-size: 1.5rem;
            color: #333;
            margin-bottom: 10px;
        }
        .stTable {
            margin-top: 1rem;
        }
    </style>
""", unsafe_allow_html=True)

# App title and description
st.markdown('<div class="header">🔬 Imprinting Factor (IF) Prediction</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">🔍 Search by Template (Compound Name) or SMILES to view IF and details.</div>', unsafe_allow_html=True)

# Input section
input_query = st.text_input(
    "Enter Template (Compound Name) or SMILES",
    placeholder="E.g., compound name or CC(=O)OC1=CC=CC=C1C(=O)O"
)

def makeblock(smi):
    """Convert SMILES to a 3D MolBlock."""
    # Split multiple SMILES if combined with '.'
    smiles_list = smi.split('.')
    
    # Combine molecules
    combined_mol = None
    for single_smile in smiles_list:
        mol = Chem.MolFromSmiles(single_smile.strip())
        if mol is None:
            raise ValueError(f"Invalid SMILES: {single_smile}")
        
        mol = Chem.AddHs(mol)
        
        # Embed the molecule
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)  # Additional optimization
        
        # Combine molecules if multiple SMILES
        if combined_mol is None:
            combined_mol = mol
        else:
            combined_mol = Chem.CombineMols(combined_mol, mol)
    
    mol_block = Chem.MolToMolBlock(combined_mol)
    return mol_block

def render_mol(xyz):
    """
    Render a molecule in 3D using py3Dmol with ball and stick representation 
    and distinct colors for Template, Monomer, and Crosslinker.
    
    Args:
        xyz (str): Molecular block or XYZ coordinates
    
    Returns:
        str: HTML representation of the 3D molecular view
    """
    viewer = py3Dmol.view(width=800, height=400)
    
    # Split the molecule into separate fragments
    mol = Chem.MolFromMolBlock(xyz)
    fragments = Chem.GetMolFrags(mol, asMols=True)
    
    # Define color and style scheme
    colors = {
        'template': 'green',     # Template in green
        'monomer': 'blue',       # Monomer in blue
        'crosslinker': 'red'     # Crosslinker in red
    }
    
    # Color each fragment with a different style and color
    for i, frag in enumerate(fragments):
        # Convert fragment to block
        frag_block = Chem.MolToMolBlock(frag)
        
        # Assign color based on fragment index
        color = colors.get(list(colors.keys())[i], 'gray')
        
        # Add fragment with ball and stick representation
        viewer.addModel(frag_block, "mol")
        viewer.setStyle({'model': -1}, {
            'stick': {
                'color': color, 
                'radius': 0.05,  # Significantly reduced stick radius
                'transparency': 0.2  # Slight transparency for better visibility
            },
            'sphere': {
                'color': color, 
                'scale': 0.3,  # Reduced sphere size
                'transparency': 0.3  # Slight transparency for atoms
            }
        })
    
    # Set background and camera settings
    viewer.setBackgroundColor("white")
    viewer.zoomTo()
    viewer.spin(True)  # Optional: Add rotation to the model
    viewer.setViewStyle({
        'style': 'outline'  # Add an outline to improve depth perception
    })

    # Return the HTML representation of the 3D view
    return viewer._make_html()

if input_query:
    # Filter dataset for matches
    matches = dataset[
        (dataset['Compound name'].str.lower() == input_query.lower()) |
        (dataset['template smile'].str.lower() == input_query.lower()) |
        (dataset['momomer smile'].str.lower() == input_query.lower()) |
        (dataset['crosslinker smile'].str.lower() == input_query.lower())
    ]

    if matches.empty:
        st.warning("⚠️ No matching Template or SMILES found.")
    else:
        # Create a dropdown to select the IF value from the matches
        if_values = matches[['Compound name', 'IF']].drop_duplicates()
        selected_if = st.selectbox("Select an IF value", if_values['IF'].unique())

        # Filter the dataset for the selected IF value
        selected_data = matches[matches['IF'] == selected_if].iloc[0]

        # Display IF Value
        st.markdown("###  Retrieved Imprinting Factor (IF) Value")
        st.markdown(f'<div class="if-value">IF Value: {selected_data["IF"]:.2f}</div>', unsafe_allow_html=True)

        # Display Template (Compound Name), Crosslinker, and Functional Monomer
        st.markdown('<div class="details-section">', unsafe_allow_html=True)
        st.markdown('<div class="details-title" style="color: #FFD700;">📋 Compound name and Related Details</div>', unsafe_allow_html=True)
        st.write(f"**🧪 Template:** {selected_data['Compound name']}")
        st.write(f"**🧬 Functional Monomer:** {selected_data['FUNCTIONAL MONOMER']}")
        st.write(f"**🔗 Crosslinker:** {selected_data['crosslinker']}")
        st.write(f"**⚗️ Template Concentration:** {selected_data['tempConc(mmol)']} mmol")
        st.write(f"**⚗️ Monomer Concentration:** {selected_data['MonomerConc(mmol)']} mmol")
        st.write(f"**⚗️ Crosslinker Concentration:** {selected_data['CROSSLINKER CONC(mmol)']} mmol")
        st.markdown('</div>', unsafe_allow_html=True)

        # Combine Template, Monomer, and Crosslinker SMILES
        combined_smiles = selected_data["template smile"] + '.' + selected_data["momomer smile"] + '.' + selected_data["crosslinker smile"]

       
        # Render 3D model for the combined SMILES
        st.markdown("### 🧪 3D Molecular Structure of Template, Monomer, and Crosslinker")
        try:
            mol_block_combined = makeblock(combined_smiles)
            mol_html_combined = render_mol(mol_block_combined)
            st.components.v1.html(mol_html_combined, height=450, scrolling=False)
        except Exception as e:
            st.warning(f"Could not generate 3D model for the combined molecule: {e}")