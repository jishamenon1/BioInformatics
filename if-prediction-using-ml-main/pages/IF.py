import streamlit as st
import pandas as pd

# Load dataset
dataset = pd.read_csv('your_dataset.csv')  # Replace with your dataset file path
dataset.columns = dataset.columns.str.strip()  # Clean column names

# Inject custom CSS for enhanced styling
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
            color: #FF6347;  /* Tomato color for IF value */
            background-color: #FFE4E1; /* Light pink background for highlight */
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
st.markdown('<div class="input-section">', unsafe_allow_html=True)
input_query = st.text_input("Enter Template (Compound Name) or SMILES", placeholder="E.g., compound name or CC(=O)OC1=CC=CC=C1C(=O)O")
st.markdown('</div>', unsafe_allow_html=True)

if input_query:
    # Filter dataset for matches
    matches = dataset[
        (dataset['Compound name'].str.lower() == input_query.lower()) |
        (dataset['template smile'].str.lower() == input_query.lower()) |
        (dataset['crosslinker smile'].str.lower() == input_query.lower()) |
        (dataset['momomer smile'].str.lower() == input_query.lower())
    ]

    if matches.empty:
        st.warning("⚠️ No matching Template or SMILES found.")
    else:
        if len(matches) > 1:
            st.write(f"⚠️ **Multiple entries ({len(matches)}) found for '{input_query}'. Please select one:**")
            # Dropdown for multiple entries
            selected_index = st.selectbox(
                "Choose an option:",
                matches.index.tolist(),
                format_func=lambda idx: f"select - IF: {matches.loc[idx, 'IF']:.2f}"
            )
            selected_data = matches.loc[selected_index]
        else:
            # Single match
            selected_data = matches.iloc[0]

        # Display IF Value
        st.markdown("###  Retrieved Imprinting Factor (IF)")
        st.markdown(f'<div class="if-value">IF Value: {selected_data["IF"]:.2f}</div>', unsafe_allow_html=True)

        # Display Template (Compound Name), Crosslinker, and Functional Monomer
        st.markdown('<div class="details-section">', unsafe_allow_html=True)
        st.markdown('<div class="details-title" style="color: #FFD700;">📋 Compound name and Related Details</div>', unsafe_allow_html=True)  # Yellow color added
        st.write(f"**🧪 Template:** {selected_data['Compound name']}")
        st.write(f"**🔗 Crosslinker:** {selected_data['crosslinker']}")
        st.write(f"**🧬 Functional Monomer:** {selected_data['FUNCTIONAL MONOMER']}")
        st.markdown('</div>', unsafe_allow_html=True)

        # Display Concentration Details
        st.markdown("### 💧 Concentration Details")
        concentration_data = {
            "Concentration Type": ["Template Concentration", "Monomer Concentration", "Crosslinker Concentration"],
            "Value (mmol)": [
                selected_data["tempConc(mmol)"],
                selected_data["MonomerConc(mmol)"],
                selected_data["CROSSLINKER CONC(mmol)"]
            ]
        }
        concentration_df = pd.DataFrame(concentration_data)
        st.table(concentration_df)
