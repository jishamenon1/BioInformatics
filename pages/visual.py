import streamlit as st
import matplotlib.pyplot as plt
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from stmol import showmol
import numpy as np
import pandas as pd
from PIL import Image
import io
import base64

st.set_page_config(layout="wide", page_title="Molecular Imprinted Polymer (MIP) Visualizer")

# Custom CSS for better styling
st.markdown("""
<style>
    .main .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
    }
    h1, h2, h3 {
        color: #1E3D59;
    }
    .stButton>button {
        background-color: #1E3D59;
        color: white;
        border-radius: 5px;
        padding: 0.5rem 1rem;
        font-weight: bold;
    }
    .stButton>button:hover {
        background-color: #2E5E8A;
        border-color: #2E5E8A;
    }
    .model-box {
        border: 1px solid #ddd;
        border-radius: 10px;
        padding: 15px;
        background-color: #f8f9fa;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    }
    
    .molecule-label {
        text-align: center;
        font-weight: bold;
        margin-bottom: 5px;
        color: #1E3D59;
    }
   
</style>
""", unsafe_allow_html=True)

st.title("Molecular Imprinted Polymer (MIP) 3D Visualizer")

# Create layout with columns for inputs
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("Input Parameters")
    # SMILES inputs
    monomer_smiles = st.text_input("Monomer SMILES", value="CC(=O)OCC(C=C)OC(=O)C", 
                                  help="Example: Methacrylic acid ester (CC(=O)OCC(C=C)OC(=O)C)")
    
    crosslinker_smiles = st.text_input("Crosslinker SMILES", value="C=C(CC(=O)OCC)C(=O)OCC", 
                                      help="Example: EGDMA (C=C(CC(=O)OCC)C(=O)OCC)")
    
    template_smiles = st.text_input("Template SMILES", value="C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)(C(=O)CO)O)N)O", 
                                   help="Example: Doxorubicin")
    
    # System parameters
    st.markdown("### System Parameters")
    monomer_concentration = st.number_input("Monomer Concentration (mmol)", value=1.0, min_value=0.1, step=0.1)
    crosslinker_concentration = st.number_input("Crosslinker Concentration (mmol)", value=25.0, min_value=0.1, step=0.1)
    template_concentration = st.number_input("Template Concentration (mmol)", value=0.46, min_value=0.01, step=0.01)
    polymerization_volume = st.number_input("Polymerization Volume (L)", value=20.4, min_value=0.1, step=0.1)
    monomer_molecular_weight = st.number_input("Monomer Molecular Weight (g/mol)", value=86.0, min_value=1.0, step=0.1)
    
    # Simulation parameters
    st.markdown("### Simulation Parameters")
    num_monomers = st.slider("Number of Monomers", min_value=5, max_value=50, value=20, step=1)
    num_crosslinkers = st.slider("Number of Crosslinkers", min_value=2, max_value=25, value=10, step=1)
    
    generate_button = st.button("Generate MIP Visualization")

# Function to process SMILES and create 3D visualization
def process_mip(monomer_smiles, crosslinker_smiles, template_smiles, 
               monomer_concentration, crosslinker_concentration, template_concentration,
               polymerization_volume, monomer_molecular_weight, num_monomers, num_crosslinkers):
    try:
        # Convert SMILES to RDKit Molecules
        monomer = Chem.MolFromSmiles(monomer_smiles)
        if monomer is None:
            return None, "Invalid monomer SMILES"
            
        crosslinker = Chem.MolFromSmiles(crosslinker_smiles)
        if crosslinker is None:
            return None, "Invalid crosslinker SMILES"
            
        template = Chem.MolFromSmiles(template_smiles)
        if template is None:
            return None, "Invalid template SMILES"
        
        # Add Hydrogens & Generate 3D Conformations
        monomer = Chem.AddHs(monomer)
        crosslinker = Chem.AddHs(crosslinker)
        template = Chem.AddHs(template)
        
        AllChem.EmbedMolecule(monomer, AllChem.ETKDG())
        AllChem.EmbedMolecule(crosslinker, AllChem.ETKDG())
        AllChem.EmbedMolecule(template, AllChem.ETKDG())
        
        # Create Polymer Matrix
        polymer = Chem.RWMol()
        
        monomer_atoms = []
        crosslinker_atoms = []
        
        # Insert Monomers into Polymer
        for i in range(num_monomers):
            start_idx = polymer.GetNumAtoms()
            polymer.InsertMol(monomer)
            monomer_atoms.extend(range(start_idx, start_idx + monomer.GetNumAtoms()))
        
        # Insert Crosslinkers into Polymer
        for i in range(num_crosslinkers):
            start_idx = polymer.GetNumAtoms()
            polymer.InsertMol(crosslinker)
            crosslinker_atoms.extend(range(start_idx, start_idx + crosslinker.GetNumAtoms()))
        
        # Insert Template into Polymer
        start_idx = polymer.GetNumAtoms()
        polymer.InsertMol(template)
        template_atoms = list(range(start_idx, start_idx + template.GetNumAtoms()))
        
        # Generate 3D Coordinates for Polymer
        AllChem.EmbedMolecule(polymer, AllChem.ETKDG())
        
        # Convert to MOL Format
        mol_block = Chem.MolToMolBlock(polymer)
        
        # Calculate Bond Lengths
        bond_lengths = []
        bond_data = []
        conf = polymer.GetConformer()
        
        for bond in polymer.GetBonds():
            idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            atom1 = polymer.GetAtomWithIdx(idx1)
            atom2 = polymer.GetAtomWithIdx(idx2)
            pos1 = np.array(conf.GetAtomPosition(idx1))
            pos2 = np.array(conf.GetAtomPosition(idx2))
            distance = np.linalg.norm(pos1 - pos2)
            bond_lengths.append(distance)
            
            component = "Monomer"
            if idx1 in crosslinker_atoms or idx2 in crosslinker_atoms:
                component = "Crosslinker"
            if idx1 in template_atoms or idx2 in template_atoms:
                component = "Template"
                
            bond_data.append({
                'Bond': f"{atom1.GetSymbol()}{idx1} - {atom2.GetSymbol()}{idx2}",
                'Length (Å)': round(distance, 2),
                'Component': component
            })
        
        # Calculate Binding Capacity - Modified as per request
        m_polymer = (monomer_concentration * num_monomers + crosslinker_concentration * num_crosslinkers) * monomer_molecular_weight
        Q_MIP = (template_concentration * polymerization_volume) / m_polymer * 1e3
        Q_NIP = template_concentration * Q_MIP
        IF = Q_MIP / Q_NIP
        
        # Create 3D Visualization
        viewer = py3Dmol.view(width=800, height=600)
        viewer.addModel(mol_block, "mol")
        
        # Ball-and-Stick Representation
        for atom_idx in monomer_atoms:
            viewer.setStyle({"model": -1, "serial": atom_idx + 1}, {"stick": {"radius": 0.25, "color": "#A7C7E7"}})
        
        for atom_idx in crosslinker_atoms:
            viewer.setStyle({"model": -1, "serial": atom_idx + 1}, {"stick": {"radius": 0.25, "color": "#F4A988"}})
        
        for atom_idx in template_atoms:
            viewer.setStyle({"model": -1, "serial": atom_idx + 1}, {"stick": {"radius": 0.3, "color": "#88D8B0"}})
        
        viewer.setBackgroundColor('white')
        viewer.zoomTo()
        
        return {
            'viewer': viewer,
            'bond_data': pd.DataFrame(bond_data),
            'bond_lengths': bond_lengths,
            'Q_MIP': Q_MIP,
            'Q_NIP': Q_NIP,
            'IF': IF,
            'polymer': polymer,
            'monomer': monomer,
            'crosslinker': crosslinker,
            'template': template
        }, None
    
    except Exception as e:
        return None, f"Error processing molecules: {str(e)}"

# When generate button is clicked
if generate_button:
    with st.spinner("Generating MIP visualization..."):
        results, error = process_mip(
            monomer_smiles, crosslinker_smiles, template_smiles,
            monomer_concentration, crosslinker_concentration, template_concentration,
            polymerization_volume, monomer_molecular_weight, num_monomers, num_crosslinkers
        )
        
        if error:
            st.error(error)
        else:
            # Main 3D visualization below input parameters
            st.subheader("3D Visualization of MIP System")
            st.markdown("#### Color code: Blue - Monomer | Orange - Crosslinker | Green - Template")
            
            # Simple viewer container
            st.markdown('<div class="viewer-container">', unsafe_allow_html=True)
            
            # Display the 3D visualization
            showmol(results['viewer'], height=600, width=800)
            
            st.markdown('</div>', unsafe_allow_html=True)
            
            # Create new column layout for remaining visualizations
            viz_col1, viz_col2 = st.columns([1, 1])
            
            with viz_col1:
                # Individual molecule visualizations
                st.subheader("Individual Molecules")
                mol_col1, mol_col2, mol_col3 = st.columns(3)
                
                # Helper function for individual molecules with enhanced styling
                def show_single_molecule(mol, title, color):
                    mol_viewer = py3Dmol.view(width=250, height=250)
                    mol_viewer.addModel(Chem.MolToMolBlock(mol), "mol")
                    mol_viewer.setStyle({}, {"stick": {"radius": 0.2, "color": color}})
                    # Add atom labels
                    for atom in mol.GetAtoms():
                        idx = atom.GetIdx()
                        mol_viewer.addLabel(atom.GetSymbol(), {"position": {"model": 0, "serial": idx + 1}, 
                                                            "backgroundColor": "transparent", 
                                                            "fontColor": "#333333", 
                                                            "fontSize": 12})
                    mol_viewer.zoomTo()
                    mol_viewer.setBackgroundColor('white')
                    return mol_viewer
                
                # Create styled containers for individual molecules
                with mol_col1:
                    st.markdown('<div class="molecule-label">Monomer</div>', unsafe_allow_html=True)
                    st.markdown('<div class="model-box">', unsafe_allow_html=True)
                    showmol(show_single_molecule(results['monomer'], "Monomer", "#A7C7E7"), height=250, width=250)
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                with mol_col2:
                    st.markdown('<div class="molecule-label">Crosslinker</div>', unsafe_allow_html=True)
                    st.markdown('<div class="model-box">', unsafe_allow_html=True)
                    showmol(show_single_molecule(results['crosslinker'], "Crosslinker", "#F4A988"), height=250, width=250)
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                with mol_col3:
                    st.markdown('<div class="molecule-label">Template</div>', unsafe_allow_html=True)
                    st.markdown('<div class="model-box">', unsafe_allow_html=True)
                    showmol(show_single_molecule(results['template'], "Template", "#88D8B0"), height=250, width=250)
                    st.markdown('</div>', unsafe_allow_html=True)

                # Display binding capacity with enhanced styling
                st.subheader("Binding Capacity Calculation")
                
                # Create a container with a background color
                st.markdown("""
                <style>
                .capacity-container {
                    background-color: #f9f9f9;
                    border-radius: 10px;
                    padding: 15px;
                    margin: 10px 0;
                    border: 1px solid #ddd;
                }
                </style>
                <div class="capacity-container">
                """, unsafe_allow_html=True)
                
                binding_col1, binding_col2, binding_col3 = st.columns(3)
                
                with binding_col1:
                    st.metric("MIP Binding Capacity", f"{results['Q_MIP']:.2f} µmol/g", 
                             delta="Reference Value", delta_color="off")
                
                with binding_col2:
                    st.metric("NIP Binding Capacity", f"{results['Q_NIP']:.2f} µmol/g", 
                             delta=f"{results['Q_NIP']/results['Q_MIP']*100:.1f}% of MIP", delta_color="off")
                
                with binding_col3:
                    st.metric("Imprinting Factor (IF)", f"{results['IF']:.2f}", 
                             delta="MIP/NIP Ratio", delta_color="off")
                
                # Formula explanation - Updated as per request
                st.markdown("""
                </div>
                <div style="font-size: 0.9em; color: #666; margin-top: 5px;">
                    <b>Formulas:</b><br>
                    • Q_MIP = (template_concentration × polymerization_volume) / (m_polymer) × 1000<br>
                    • Q_NIP = template_concentration × Q_MIP<br>
                    • IF = Q_MIP / Q_NIP
                </div>
                """, unsafe_allow_html=True)
                
            with viz_col2:
                # Bond length histogram with enhanced visualization
                st.subheader("Bond Length Distribution")
                
                # Create figure with better styling
                fig, ax = plt.subplots(figsize=(8, 4))
                
                # Create separate histograms for different components
                monomer_bonds = [row['Length (Å)'] for i, row in results['bond_data'].iterrows() if row['Component'] == 'Monomer']
                crosslinker_bonds = [row['Length (Å)'] for i, row in results['bond_data'].iterrows() if row['Component'] == 'Crosslinker']
                template_bonds = [row['Length (Å)'] for i, row in results['bond_data'].iterrows() if row['Component'] == 'Template']
                
                # Plot histograms with transparency for overlapping visibility
                if monomer_bonds:
                    ax.hist(monomer_bonds, bins=15, color="#A7C7E7", edgecolor="black", alpha=0.6, label='Monomer')
                if crosslinker_bonds:
                    ax.hist(crosslinker_bonds, bins=15, color="#F4A988", edgecolor="black", alpha=0.6, label='Crosslinker')
                if template_bonds:
                    ax.hist(template_bonds, bins=15, color="#88D8B0", edgecolor="black", alpha=0.6, label='Template')
                
                # Add statistics to the plot
                all_bonds = results['bond_lengths']
                ax.axvline(np.mean(all_bonds), color='red', linestyle='dashed', linewidth=1, label=f'Mean: {np.mean(all_bonds):.2f} Å')
                
                # Add styling
                ax.set_xlabel("Bond Length (Å)", fontsize=12)
                ax.set_ylabel("Frequency", fontsize=12)
                ax.set_title("Bond Length Distribution by Component", fontsize=14)
                ax.grid(True, linestyle='--', alpha=0.7)
                ax.legend()
                
                # Set background color
                fig.patch.set_facecolor('#f8f9fa')
                ax.set_facecolor('#f8f9fa')
                
                st.pyplot(fig)
                
                # Display statistics in a box
                st.markdown(f"""
                <div style="background-color: #f0f0f0; padding: 10px; border-radius: 5px; margin-top: 10px;">
                    <b>Bond Length Statistics:</b><br>
                    • Mean: {np.mean(all_bonds):.2f} Å<br>
                    • Minimum: {np.min(all_bonds):.2f} Å<br>
                    • Maximum: {np.max(all_bonds):.2f} Å<br>
                    • Standard Deviation: {np.std(all_bonds):.2f} Å
                </div>
                """, unsafe_allow_html=True)
                
            # Bond length table with search and filter capabilities (full width)
            st.subheader("Bond Length Details")
            
            # Add filter options
            filter_col1, filter_col2 = st.columns([1, 2])
            
            with filter_col1:
                component_filter = st.selectbox(
                    "Filter by Component",
                    ["All Components", "Monomer", "Crosslinker", "Template"]
                )
            
            with filter_col2:
                search_term = st.text_input("Search Bonds (e.g., 'C-O', 'N')", "")
            
            # Apply filters
            filtered_data = results['bond_data'].copy()
            
            if component_filter != "All Components":
                filtered_data = filtered_data[filtered_data['Component'] == component_filter]
            
            if search_term:
                filtered_data = filtered_data[filtered_data['Bond'].str.contains(search_term, case=False)]
            
            # Display dataframe with enhanced styling
            st.dataframe(
                filtered_data,
                column_config={
                    "Bond": st.column_config.TextColumn(
                        "Bond",
                        help="Atom pair forming the bond"
                    ),
                    "Length (Å)": st.column_config.NumberColumn(
                        "Length (Å)",
                        format="%.2f Å",
                        help="Distance between atoms in Angstroms"
                    ),
                    "Component": st.column_config.SelectboxColumn(
                        "Component",
                        options=["Monomer", "Crosslinker", "Template"],
                        width="medium",
                        help="Molecular component where this bond is found"
                    )
                },
                hide_index=True,
                use_container_width=True
            )
            
            # Add download button for bond data
            csv = filtered_data.to_csv(index=False)
            b64 = base64.b64encode(csv.encode()).decode()
            href = f'<a href="data:file/csv;base64,{b64}" download="bond_length_data.csv" style="display: inline-block; padding: 0.5rem 1rem; background-color: #1E3D59; color: white; text-decoration: none; border-radius: 4px; margin-top: 10px;">Download Bond Data as CSV</a>'
            st.markdown(href, unsafe_allow_html=True)

# Add information about the application
st.sidebar.title("About")

st.sidebar.subheader("Example Molecules")
st.sidebar.markdown("""
**Monomers:**
- Methacrylic acid ester: `CC(=O)OCC(C=C)OC(=O)C`
- Acrylic acid: `C=CC(=O)O`

**Crosslinkers:**
- EGDMA: `C=C(CC(=O)OCC)C(=O)OCC`
- TRIM: `CC(=C)C(=O)OCC(COC(=O)C(=C)C)OC(=O)C(=C)C`

**Templates:**
- Doxorubicin: `C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)(C(=O)CO)O)N)O`
- Caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
""")

st.sidebar.subheader("Legend")
st.sidebar.markdown("""
- **Blue**: Monomer molecules
- **Orange**: Crosslinker molecules
- **Green**: Template molecule
""")