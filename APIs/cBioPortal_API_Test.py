import requests

BASE = "https://www.cbioportal.org/api"
GENE = "PTDSS1"
STUDY = "paad_tcga_pan_can_atlas_2018"


# Send a request and return the answer as Python data.
def get(path, **params):
    response = requests.get(BASE + path, params=params, timeout=60)
    response.raise_for_status()
    return response.json()


# 1. Look up a gene.
print("\nGENE INFORMATION")
print(get(f"/genes/{GENE}"))

# 2. Get information about the pancreatic cancer study.
print("\nSTUDY INFORMATION")
study = get(f"/studies/{STUDY}")
print("Name:", study["name"])
print("Description:", study.get("description", "Not provided"))

# 3. See which measurements are available.
print("\nAVAILABLE MOLECULAR PROFILES")
profiles = get(f"/studies/{STUDY}/molecular-profiles")
for profile in profiles:
    print(profile["molecularProfileId"], "-", profile["name"])

# 4. List sample groups, such as sequenced samples.
print("\nSAMPLE GROUPS")
sample_lists = get(f"/studies/{STUDY}/sample-lists")
for sample_list in sample_lists:
    print(sample_list["sampleListId"], "-", sample_list["name"])

# 5. Count sequenced samples (not necessarily unique patients).
print("\nSEQUENCED SAMPLE COUNT")
sample_ids = get(f"/sample-lists/{STUDY}_sequenced/sample-ids")
print("Sequenced samples:", len(sample_ids))

# 6. List the clinical fields available in this study.
print("\nAVAILABLE CLINICAL INFORMATION")
attributes = get(f"/studies/{STUDY}/clinical-attributes")
for attribute in attributes:
    print(attribute["clinicalAttributeId"], "-", attribute["displayName"])
