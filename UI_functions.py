import streamlit as st
from Standard_Accretion_Disk import get_region_a_data
from The_Plotting_function import plot_region_a
from constants import R_AB_LOOKUP

# USER DEFINED FUNCTIONS FOR USER INTERFACE 

def display_img(image_file, preview_width=400):

    ext = os.path.splitext(image_file)[1].lower()

    if ext in [".jpg", ".jpeg"]:
        mime = "jpeg"
    elif ext == ".png":
        mime = "png"
    else:
        mime = "jpeg"   # safe fallback

    with open(image_file, "rb") as f:
        data = f.read()

    encoded = base64.b64encode(data).decode()

def render_region_a_ui():
    """Streamlit layout component for Region (a)."""
    st.sidebar.subheader("Region (a) Parameters")
    m_val = st.sidebar.selectbox("Black Hole Mass ($m$)", [1e8, 5e8, 1e9], index=0)
    mdot_val = st.sidebar.selectbox("Accretion Rate ($\dot{m}$)", [0.01, 0.05], index=0)
    alpha_val = st.sidebar.number_input("Alpha ($\alpha$)", value=0.1, step=0.01)
    f1_val = st.sidebar.number_input("$f_1$ factor", value=1.0, step=0.1)
    eta_E_val = st.sidebar.number_input("Efficiency ($\eta_E$)", value=0.06, step=0.01)

    # Compute and plot
    r_grid, T_vals, n_vals, r_ab = get_region_a_data(
        m_val, mdot_val, alpha_val, f1_val, eta_E_val
    )
    
    st.info(f"Calculated Boundary Radius $r_{{s,ab}} = {r_ab} \, R_S$")
    fig = plot_region_a(r_grid, T_vals, n_vals, r_ab)
    st.pyplot(fig)

def render_sidebar():
    # Existing M and mdot inputs...
    m_val = st.sidebar.selectbox("Black Hole Mass (M/M_sun)", [1e8, 5e8, 1e9])
    mdot_val = st.sidebar.selectbox("Accretion Rate (m_dot)", [0.01, 0.05])
    
    # New restored parameters
    alpha_val = st.sidebar.number_input("Viscosity (α)", value=0.1)
    f1_val = st.sidebar.number_input("f1 factor", value=1.0)
    eta_E_val = st.sidebar.number_input("Efficiency (η_E)", value=0.06)
    
    r_ab = R_AB_LOOKUP.get((m_val, mdot_val), 50.54)
    st.sidebar.info(f"r_ab boundary: {r_ab} R_S")
    
    return m_val, mdot_val, alpha_val, f1_val, eta_E_val, r_ab
    
    st.markdown(
        f"""
        <div style="text-align:center;">
            <img src="data:image/{mime};base64,{encoded}" 
                 style="max-width:90%; width:{preview_width}px; 
                        border:4px solid #ccc; 
                        box-shadow:5px 5px 15px rgba(0,0,0,0.3); 
                        border-radius:8px;">
        </div>
        """,
        unsafe_allow_html=True
    )
