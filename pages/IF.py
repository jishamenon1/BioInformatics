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

def makeblock(smi):
    """Convert SMILES to a 3D MolBlock with improved error handling."""
    try:
        # Split multiple SMILES if combined with '.'
        smiles_list = smi.split('.')
        
        # Combine molecules
        combined_mol = None
        for single_smile in smiles_list:
            mol = Chem.MolFromSmiles(single_smile.strip())
            if mol is None:
                raise ValueError(f"Invalid SMILES: {single_smile}")
            
            mol = Chem.AddHs(mol)
            
            try:
                # Try embedding with different random seeds if initial embedding fails
                for seed in [42, 123, 987, 555]:
                    if AllChem.EmbedMolecule(mol, randomSeed=seed) == 0:
                        AllChem.MMFFOptimizeMolecule(mol)
                        break
                else:
                    raise ValueError(f"Could not embed molecule: {single_smile}")
            except Exception as e:
                raise ValueError(f"Error in 3D conformation generation: {str(e)}")
            
            # Combine molecules if multiple SMILES
            if combined_mol is None:
                combined_mol = mol
            else:
                combined_mol = Chem.CombineMols(combined_mol, mol)
        
        if combined_mol is None:
            raise ValueError("No valid molecules to process")
            
        mol_block = Chem.MolToMolBlock(combined_mol)
        return mol_block
        
    except Exception as e:
        raise ValueError(f"Error generating molecular block: {str(e)}")

def render_mol(xyz):
    """
    Render a molecule in 3D using py3Dmol with ball and stick representation 
    and distinct colors for Template, Monomer, and Crosslinker.
    """
    viewer = py3Dmol.view(width=800, height=400)
    
    try:
        # Split the molecule into separate fragments
        mol = Chem.MolFromMolBlock(xyz)
        if mol is None:
            raise ValueError("Invalid molecular block")
            
        fragments = Chem.GetMolFrags(mol, asMols=True)
        
        # Define color scheme with fallback colors
        colors = [
            '#2ECC71',  # Green for template
            '#3498DB',  # Blue for monomer
            '#E74C3C',  # Red for crosslinker
            '#F1C40F',  # Yellow (fallback)
            '#9B59B6',  # Purple (fallback)
            '#1ABC9C',  # Turquoise (fallback)
        ]
        
        # Add each fragment with its own color
        for i, frag in enumerate(fragments):
            # Convert fragment to block
            frag_block = Chem.MolToMolBlock(frag)
            
            # Use modulo to cycle through colors if there are more fragments than colors
            color = colors[i % len(colors)]
            
            # Add fragment with ball and stick representation
            viewer.addModel(frag_block, "mol")
            viewer.setStyle({'model': -1}, {
                'stick': {
                    'color': color, 
                    'radius': 0.05,
                    'transparency': 0.2
                },
                'sphere': {
                    'color': color, 
                    'scale': 0.3,
                    'transparency': 0.3
                }
            })
        
        # Set background and camera settings
        viewer.setBackgroundColor("white")
        viewer.zoomTo()
        viewer.spin(True)
        viewer.setViewStyle({'style': 'outline'})
        
        return viewer._make_html()
        
    except Exception as e:
        raise ValueError(f"Error rendering 3D model: {str(e)}")

# App title and description
st.markdown('<div class="header">🔬 Imprinting Factor (IF) Prediction</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">🔍 Search by Template (Compound Name) or SMILES to view IF and details.</div>', unsafe_allow_html=True)

# Input section
input_query = st.text_input(
    "Enter Template (Compound Name) or SMILES",
    placeholder="E.g., compound name or CC(=O)OC1=CC=CC=C1C(=O)O"
)

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
        combined_smiles = f"{selected_data['template smile']}.{selected_data['momomer smile']}.{selected_data['crosslinker smile']}"

        # Display 3D model with improved error handling
        st.markdown("### 🧪 3D Molecular Structure of Template, Monomer, and Crosslinker")
        
        try:
            mol_block_combined = makeblock(combined_smiles)
            mol_html_combined = render_mol(mol_block_combined)
            st.components.v1.html(mol_html_combined, height=450, scrolling=False)
        except Exception as e:
            st.error(f"⚠️ 3D Model Generation Error: {str(e)}")
            st.info("💡 This might happen due to complex molecular structures or invalid SMILES strings. The rest of the data is still valid.")
            
            # Display the SMILES strings for debugging
            with st.expander("🔍 Show SMILES strings for debugging"):
                st.write("Template SMILES:", selected_data['template smile'])
                st.write("Monomer SMILES:", selected_data['momomer smile'])
                st.write("Crosslinker SMILES:", selected_data['crosslinker smile'])
