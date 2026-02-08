import requests
import time

class ACMGClassifier:
    def __init__(self):
        self.pm2_threshold = 0.0001
        self.lof_consequences = ['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant']

    def get_variant_data(self, variant_id):
        """Ensembl API - Biological Prediction"""
        server = 'https://rest.ensembl.org'
        ext = f"/vep/human/hgvs/{variant_id}?"
        try:
            response = requests.get(server + ext, headers={"Content-Type": "application/json"})
            return response.json()[0] if response.ok else None
        except: return None

    def get_clinvar_data(self, variant_id):
        """NCBI API - Human Clinical History"""
        # Step 1: Search for the ID
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={variant_id}&retmode=json"
        try:
            search_res = requests.get(search_url).json()
            id_list = search_res.get('esearchresult', {}).get('idlist', [])
            
            if not id_list:
                return "No ClinVar Record Found"

            # Step 2: Get the summary for that ID
            summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={id_list[0]}&retmode=json"
            summary_res = requests.get(summary_url).json()
            return summary_res['result'][id_list[0]]['clinical_significance']['description']
        except:
            return "ClinVar Search Failed"

    def calculate_logic(self, data):
        evidence = {'pvs': 0, 'pm': 0}
        if not data: return evidence, "unknown", 0.0
        
        colocated = data.get('colocated_variants', [{}])[0]
        freq = colocated.get('gnomad_af', 0.0)
        if freq < self.pm2_threshold: evidence['pm'] += 1

        consequence = data.get('most_severe_consequence', '')
        if consequence in self.lof_consequences: evidence['pvs'] += 1
            
        return evidence, consequence, freq

    def get_final_verdict(self, scores, clinvar_status):
        pvs, pm = scores['pvs'], scores['pm']
        
        # Calculate our internal prediction
        if pvs >= 1 and pm >= 1: prediction = "PATHOGENIC"
        elif pvs >= 1: prediction = "LIKELY PATHOGENIC"
        else: prediction = "VUS"

        # Check for conflicts
        if prediction == "PATHOGENIC" and "Pathogenic" in clinvar_status:
            return f"CONFIRMED {prediction}"
        elif prediction == "PATHOGENIC" and "Benign" in clinvar_status:
            return "WARNING: CONFLICTING DATA (Predict: Path, ClinVar: Benign)"
        
        return f"{prediction} (ClinVar: {clinvar_status})"

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    app = ACMGClassifier()
    variant = "NM_007294.4:c.5123C>G" 
    
    # Fetch from two different sources
    ens_data = app.get_variant_data(variant)
    clinvar_report = app.get_clinvar_data(variant)
    
    # Process
    results, cons, freq = app.calculate_logic(ens_data)
    verdict = app.get_final_verdict(results, clinvar_report)
    
    print("-" * 50)
    print(f"ANALYSIS FOR:  {variant}")
    print(f"CONSEQUENCE:   {cons}")
    print(f"CLINVAR SAYS:  {clinvar_report}")
    print(f"VERDICT:       {verdict}")
    print("-" * 50)


            


    
    
    
    
    
