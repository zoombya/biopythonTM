import streamlit as st
import pandas as pd
import altair as alt
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

# Optimize screen usage
st.set_page_config(layout="wide")

# --- Zone 1: Header (Input file and Method Selection) ---
with st.container():
    st.title("DNA Melting Temperature Calculator")
    st.header("Input and Method Selection")
    
    # Sequence input: either manual text or file upload.
    seq_text = st.text_area("Enter nucleotide sequences (one per line):", height=100)
    uploaded_file = st.file_uploader("Or upload a text file", type=["txt"])
    if uploaded_file:
        seq_text = uploaded_file.getvalue().decode("utf-8")
    sequences = [s.strip() for s in seq_text.splitlines() if s.strip()]
    
    # Method selection (common to all methods)
    method = st.selectbox("Select Calculation Method", 
                          options=["Tm_Wallace", "Tm_GC", "Tm_NN"])

# --- Zone 2: Two Columns ---
col_output, col_params = st.columns([2, 1])
results = []
if sequences:
    # --- Parameters in Right Column ---
    with col_params:
        st.header("Parameters")
        if method in ["Tm_Wallace", "Tm_GC", "Tm_NN"]:
            check = st.checkbox("Validate sequence", value=True)
            strict = st.checkbox("Strict mode", value=True)
        
        if method == "Tm_Wallace":
            st.info("Calculates Tm = 4°C×(G+C) + 2°C×(A+T).")
        
        elif method == "Tm_GC":
            st.info("Empirical GC-based method with salt and mismatch corrections.")
            valueset = st.number_input("Valueset", min_value=1, max_value=10, value=7, step=1)
            userset_str = st.text_input("User-defined set (A,B,C,D) [optional]", value="")
            try:
                userset = tuple(float(x.strip()) for x in userset_str.split(",")) if userset_str else None
                if userset and len(userset) != 4:
                    st.error("User-defined set must contain 4 numbers.")
                    userset = None
            except Exception:
                st.error("Error parsing user-defined set.")
                userset = None
            na_conc = st.number_input("Na⁺ (mM)", min_value=0.0, value=50.0, step=0.1)
            k_conc = st.number_input("K⁺ (mM)", min_value=0.0, value=0.0, step=0.1)
            tris_conc = st.number_input("Tris (mM)", min_value=0.0, value=0.0, step=0.1)
            mg_conc = st.number_input("Mg²⁺ (mM)", min_value=0.0, value=0.0, step=0.1)
            dntp_conc = st.number_input("dNTPs (mM)", min_value=0.0, value=0.0, step=0.1)
            saltcorr = st.number_input("Salt Correction Type (0=none, 1–7)", min_value=0, max_value=7, value=5, step=1)
            mismatch = st.checkbox("Apply mismatch correction", value=True)
        
        elif method == "Tm_NN":
            st.info("Nearest Neighbor thermodynamics calculation.")
            c_seq = st.text_input("Complementary sequence (optional)", value="")
            shift = st.number_input("Shift (offset)", value=0, step=1)
            nn_table = st.selectbox("Nearest Neighbor Table", 
                                    options=["DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4",
                                             "RNA_NN1", "RNA_NN2", "RNA_NN3", "R_DNA_NN1"])
            tmm_table = st.selectbox("Terminal Mismatch Table", options=["DNA_TMM1"])
            imm_table = st.selectbox("Internal Mismatch Table", options=["DNA_IMM1"])
            de_table = st.selectbox("Dangling End Table", options=["DNA_DE1", "RNA_DE1"])
            dnac1 = st.number_input("DNA Concentration 1 (nM)", min_value=0.0, value=25.0, step=0.1)
            dnac2 = st.number_input("DNA Concentration 2 (nM)", min_value=0.0, value=25.0, step=0.1)
            selfcomp = st.checkbox("Self-complementary", value=False)
            na_conc = st.number_input("Na⁺ (mM)", min_value=0.0, value=50.0, step=0.1)
            k_conc = st.number_input("K⁺ (mM)", min_value=0.0, value=0.0, step=0.1)
            tris_conc = st.number_input("Tris (mM)", min_value=0.0, value=0.0, step=0.1)
            mg_conc = st.number_input("Mg²⁺ (mM)", min_value=0.0, value=0.0, step=0.1)
            dntp_conc = st.number_input("dNTPs (mM)", min_value=0.0, value=0.0, step=0.1)
            saltcorr = st.number_input("Salt Correction Type (0=none, 1–7)", min_value=0, max_value=7, value=5, step=1)
            nn_table_obj = getattr(mt, nn_table)
            tmm_table_obj = getattr(mt, tmm_table)
            imm_table_obj = getattr(mt, imm_table)
            de_table_obj = getattr(mt, de_table)
    
    # --- Output in Left Column ---
    with col_output:
        st.header("Output")
        for s in sequences:
            try:
                seq_obj = Seq(s)
                if method == "Tm_Wallace":
                    tm = mt.Tm_Wallace(seq_obj, check=check, strict=strict)
                elif method == "Tm_GC":
                    tm = mt.Tm_GC(
                        seq_obj,
                        check=check,
                        strict=strict,
                        valueset=valueset,
                        userset=userset,
                        Na=na_conc * 1e-3,     # Convert mM to M
                        K=k_conc * 1e-3,
                        Tris=tris_conc * 1e-3,
                        Mg=mg_conc * 1e-3,
                        dNTPs=dntp_conc * 1e-3,
                        saltcorr=saltcorr,
                        mismatch=mismatch
                    )
                elif method == "Tm_NN":
                    tm = mt.Tm_NN(
                        seq_obj,
                        check=check,
                        strict=strict,
                        c_seq=c_seq.strip() if c_seq.strip() else None,
                        shift=shift,
                        nn_table=nn_table_obj,
                        tmm_table=tmm_table_obj,
                        imm_table=imm_table_obj,
                        de_table=de_table_obj,
                        dnac1=dnac1 * 1e-9,    # Convert nM to M
                        dnac2=dnac2 * 1e-9,
                        selfcomp=selfcomp,
                        Na=na_conc * 1e-3,
                        K=k_conc * 1e-3,
                        Tris=tris_conc * 1e-3,
                        Mg=mg_conc * 1e-3,
                        dNTPs=dntp_conc * 1e-3,
                        saltcorr=saltcorr
                    )
                results.append((s, tm))
            except Exception as e:
                results.append((s, f"Error: {e}"))
        
        # Generate enlarged histogram using Altair with bins 5°C apart.
        num_tms = [tm for _, tm in results if isinstance(tm, (int, float))]
        if len(num_tms) > 1:
            df_hist = pd.DataFrame({"Tm (°C)": num_tms})
            hist_chart = alt.Chart(df_hist).mark_bar().encode(
                alt.X("Tm (°C):Q", bin=alt.Bin(step=5), title="Tm (°C)"),
                alt.Y("count()", title="Frequency")
            ).properties(title="Tm Distribution", width=600, height=400)
            st.altair_chart(hist_chart, use_container_width=True)
        
        # --- Zone 3: Full-Width Data Frame ---
        if sequences:
            df = pd.DataFrame(results, columns=["Sequence", "Tm (°C)"])
            df.sort_values("Tm (°C)", ascending=False, inplace=True)
            st.subheader("Detailed Results")
            st.dataframe(df)