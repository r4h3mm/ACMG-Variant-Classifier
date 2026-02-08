🧬 ACMG Variant Classifier
An automated bioinformatics tool designed to streamline the clinical interpretation of human genetic variants according to ACMG/AMP guidelines. This application aggregates data from global genomic databases to provide a preliminary pathogenicity verdict.

🚀 Key Features
Universal Gene Support: While initially tested on BRCA1, the engine is gene-agnostic and retrieves data for any human transcript via Ensembl APIs.

Real-time Data Integration: * Ensembl VEP: Fetches consequence (e.g., Frameshift, Missense) and population allele frequencies.

NCBI ClinVar: Cross-references variants against the world's largest database of clinical interpretations.

Automated Evidence Scoring:

PVS1 (Very Strong): Automatically triggered for null variants (frameshifts) in genes where loss-of-function is a known disease mechanism.

PM2 (Moderate): Triggered when the variant is absent or extremely rare in global populations (gnomAD).

💻 Tech Stack
Language: Python 3.9+

Frontend: Streamlit (Web Dashboard)

APIs: Ensembl REST API, NCBI E-utilities

📥 Installation & Local Setup
To run this project on your local machine:

Clone the repository:

Bash
git clone https://github.com/r4h3mm/ACMG-Variant-Classifier.git
cd ACMG-Variant-Classifier
Create and activate a virtual environment:

Bash
python3 -m venv acmg_env
source acmg_env/bin/activate
Install dependencies:

Bash
pip install -r requirements.txt
Launch the app:

Bash
streamlit run app.py
⚠️ Disclaimer
This tool is for educational and research purposes only. It should not be used for clinical diagnostics or medical decision-making. Always refer to official ClinGen or professional genetic counselor interpretations.
