import streamlit as st
import requests

class ACMGClassifier:
    def __init__(self):
        self.pm2_threshold = 0.0001
        self.lof_consequences = ['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant']

    def get_variant_data(self, variant_id):
        server = 'https://rest.ensembl.org'
        ext = f"/vep/human/hgvs/{variant_id}?"
        try:
            response = requests.get(server + ext, headers={"Content-Type": "application/json"})
            return response.json()[0] if response.ok else None
        except: return None

    def get_clinvar_data(self, variant_id):
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={variant_id}&retmode=json"
        try:
            search_res = requests.get(search_url).json()
            id_list = search_res.get('esearchresult', {}).get('idlist', [])
            if not id_list: return "No Record Found"
            summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={id_list[0]}&retmode=json"
            summary_res = requests.get(summary_url).json()
            return summary_res['result'][id_list[0]]['clinical_significance']['description']
        except: return "Search Timeout"

    def calculate_logic(self, data):
        evidence = {'pvs': 0, 'pm': 0}
        if not data: return evidence, "unknown", 0.0
        colocated = data.get('colocated_variants', [{}])[0]
        freq = colocated.get('gnomad_af', 0.0)
        if freq < self.pm2_threshold: evidence['pm'] += 1
        consequence = data.get('most_severe_consequence', '')
        if consequence in self.lof_consequences: evidence['pvs'] += 1
        return evidence, consequence, freq

# --- STREAMLIT UI ---
st.set_page_config(page_title="Genomic Variant Interpreter", page_icon="🧬")

st.title("🧬 ACMG Variant Classifier")
st.markdown("Automated ACMG Evidence Aggregator for Human Genomic Variants")

# Sidebar for Disclaimer
with st.sidebar:
    st.warning("⚠️ **Educational Tool Only**: Not for clinical use or diagnostic decisions.")
    st.info("Uses Ensembl VEP and NCBI ClinVar APIs.")

# User Input
variant_input = st.text_input("Enter Variant (HGVS format):", value="NM_007294.4:c.5137_5138insA")

if st.button("Analyze Variant"):
    app = ACMGClassifier()
    
    with st.spinner('Gathering Genomic Evidence...'):
        ens_data = app.get_variant_data(variant_input)
        clinvar_report = app.get_clinvar_data(variant_input)
        results, cons, freq = app.calculate_logic(ens_data)

    # Display Results in "Cards"
# UI Results
        st.divider()
        
        # We give the first and third columns more width (ratio 2:1:2)
        col1, col2, col3 = st.columns([2, 1, 2]) 
        
        with col1:
            st.markdown(f"**Consequence**\n\n{cons.replace('_', ' ').title()}")
            
        with col2:
            st.markdown(f"**Frequency**\n\n{freq:.6f}")
            
        with col3:
            st.markdown(f"**ClinVar Status**\n\n{clinvar_report}") 


    # Final Verdict Logic
    pvs, pm = results['pvs'], results['pm']
    if pvs >= 1 and pm >= 1:
        st.error(f"### FINAL VERDICT: PATHOGENIC")
        st.write("Criteria met: PVS1 (Very Strong) and PM2 (Moderate)")
    elif pvs >= 1:
        st.warning(f"### FINAL VERDICT: LIKELY PATHOGENIC")
    else:
        st.success(f"### FINAL VERDICT: VUS / BENIGN")
