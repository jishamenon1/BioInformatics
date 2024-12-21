import streamlit as st

# Display the main title image at the top and center it using Streamlit
st.markdown("""
    <style>
        .center-img {
            display: flex;
            justify-content: center;
        }
        .tag {
            text-align: left;
            font-size: 18px;
            color: #B60E54;  /* Set color of the tag to match the heading color */
            font-weight: bold;
            font-style: italic;
            margin-top: 10px;
        }
    </style>
""", unsafe_allow_html=True)

# Add the image using st.image() directly (no concatenation with strings)
st.image("assets/amlogo.png", width=350)

# Add the "Computational BioSciences" tag below the logo
st.markdown('<div class="tag">[TAG:- Computational BioSciences] </div>', unsafe_allow_html=True)

# Custom CSS to increase title font size, style email icons, and add color
st.markdown("""
<style>
.section-heading {
    font-size: 22px !important;  /* Increased font size for the title */
    font-weight: bold;
    margin-bottom: 20px;
    color: #B60E54;  /* Set color of the title to a shade of pink */
}

.email-icon {
    margin-right: 8px;
    color: #0066cc;  /* Set color of the email icon to blue */
}
</style>
""", unsafe_allow_html=True)

# Create two columns for left and right alignment (50%-50%)
col1, col2 = st.columns(2)

# Left column content: Amrita School of Computing
with col1:
    st.markdown("""<h3 class="section-heading">Amrita School of Computing, Amrita Vishwa Vidyapeetham, Amritapuri Campus</h3>""", unsafe_allow_html=True)
    st.markdown("""  
    **Dr. Manjusha Nair**  
    📧 Email: [manjushanair@am.amrita.edu](mailto:manjushanair@am.amrita.edu)
    
    **Ms. Jisha R C**  
    📧 Email: [Jisha@am.amrita.edu](mailto:Jisha@am.amrita.edu)
    
    **Mr. Karthik P**  
    📧 Email: [karthikadi2001@gmail.com](mailto:karthikadi2001@gmail.com)
    
    **Mr. Achu P**  
    📧 Email: [achup030801@gmail.com](mailto:achup030801@gmail.com)
    """, unsafe_allow_html=True)

# Right column content: Amrita School of Biotechnology
with col2:
    st.markdown("""<h3 class="section-heading">Amrita School of Biotechnology, Amrita Vishwa Vidyapeetham, Amritapuri Campus</h3>""", unsafe_allow_html=True)
    st.markdown("""  
    **Dr. Renuka Suravajhala**  
    📧 Email: [renus@am.amrita.edu](mailto:renus@am.amrita.edu)

    **Aswathi P**  
    📧 Email: [aswathipnambiar99@gmail.com](mailto:aswathipnambiar99@gmail.com)
    """, unsafe_allow_html=True)
