import streamlit as st
import pandas as pd

# Load dataset
dataset = pd.read_csv('your_dataset.csv')  # Replace with your dataset file path
dataset.columns = dataset.columns.str.strip()  # Clean column names

# Rename and select columns
columns_to_display = {
    "Compound name": "Template",
    "template smile": "Template SMILES",
    "FUNCTIONAL MONOMER": "Functional Monomer",
    "momomer smile": "Monomer SMILES",
    "crosslinker": "Crosslinker",
    "crosslinker smile": "Crosslinker SMILES"
}
table_data = dataset[list(columns_to_display.keys())].rename(columns=columns_to_display)

# Remove rows where any column is missing
table_data = table_data.dropna(subset=columns_to_display.values())

# Group by Template and keep the first unique row for each compound name
table_data = table_data.drop_duplicates(subset=["Template"])

# Inject custom CSS for styling
st.markdown("""
    <style>
        body {
            background-color: #1e1e1e; /* Dark background */
            color: #f4f4f4;
            font-family: 'Arial', sans-serif;
        }
        .header {
            font-size: 2rem;
            font-weight: bold;
            text-align: center;
            color: #00FA9A; /* Medium spring green */
            margin-bottom: 1rem;
        }
        .sub-header {
            text-align: center;
            color: #87CEFA; /* Light sky blue */
            font-size: 1rem;
            margin-bottom: 2rem;
        }
        .stTable {
            margin-top: -2rem;
            border: 1px solid #444;
            border-radius: 10px;
            overflow-x: auto; /* Enable horizontal scrolling */
            overflow-y: auto; /* Enable vertical scrolling */
            max-height: 500px; /* Reduced max-height */
            width: 135%; /* Reduced width */
            display: block;
            margin-left: -7rem; /* Adjusted left alignment */
        }
        .stTable th {
            background-color: #333;
            color: #FFD700; /* Gold for header text */
            font-weight: bold;
            font-size: 0.9rem; /* Reduced font size */
            padding: 4px; /* Reduced padding */
        }
        .stTable td {   
            text-align: center;
            font-size: 0.9rem; /* Reduced font size */
            color: #f4f4f4; /* Light text for readability */
            padding: 3px; /* Reduced padding */
        }
        .stTable tr:nth-child(even) {
            background-color: #2e2e2e; /* Slightly lighter dark background for alternating rows */
        }
        .stTable tr:nth-child(odd) {
            background-color: #1e1e1e;
        }
        .stTable td:first-child {
            color: #FF4500; /* OrangeRed for Template */
        }
        .stTable td:nth-child(2) {
            color: #00CED1; /* Dark Turquoise for Template SMILES */
        }
        .stTable td:nth-child(3) {
            color: #32CD32; /* Lime Green for Functional Monomer */
        }
        .stTable td:nth-child(4) {
            color: #FFA500; /* Orange for Crosslinker */
        }
        .stTable td:nth-child(5) {
            color: #8A2BE2; /* BlueViolet for Crosslinker SMILES */
        }
        .stTable td:nth-child(6) {
            color: #DC143C; /* Crimson for Monomer SMILES */
        }
    </style>
""", unsafe_allow_html=True)

# Display the table with reduced size
if table_data.empty:
    st.warning("⚠️ No valid rows to display. Ensure that all columns have values.")
else:
    # Display table with scrollable area and reduced size
    st.markdown('<div style="overflow-x:auto; overflow-y:auto; width:100%; max-width:100%; margin-left: 0; margin-top: 1rem;">', unsafe_allow_html=True)
    st.table(table_data)
    st.markdown('</div>', unsafe_allow_html=True)
