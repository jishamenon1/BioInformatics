import streamlit as st
import base64

# Set the custom page title and favicon (This must be at the very top)
st.set_page_config(
    page_title="🔬 IF Prediction",  # Custom title without "Streamlit"
    page_icon="assets/Logov.png",  # Path to the favicon image
)

# Function to convert image to base64
def image_to_base64(image_path):
    with open(image_path, "rb") as img_file:
        img_base64 = base64.b64encode(img_file.read()).decode('utf-8')
    return img_base64

# Load the background image
image_path = 'assets/home_bg.png'  # Replace with the path to your image

try:
    image_base64 = image_to_base64(image_path)

    # Embed the background image and custom styles
    page_bg_img = f"""
    <style>
    .stApp {{
      background-color: #1E1E1E; 
      background-image: url("data:image/png;base64,{image_base64}"); 
      background-position: right;  /* Adjusts the position */
      background-size: 82%;  /* Zooms out the image */
      background-repeat: no-repeat;
      padding: 3rem;
      border-radius: 10px;
      text-align: center;  
      padding-left: 80px;  /* Add space for sidebar */
    }}

    /* Adjust the sidebar width (optional) */
    .css-1g1z59i {{
      width: 200px;  /* Increase the width of the sidebar */
    }}

    /* Adjust the main content area */
    .css-1y4v5i9 {{
      margin-left: 220px;  /* Shift main content to the right */
    }}
    </style>
    """
    st.markdown(page_bg_img, unsafe_allow_html=True)
except FileNotFoundError:
    st.error(f"Image not found: {image_path}")