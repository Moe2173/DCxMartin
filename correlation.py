import requests

BASE = 'https://www.cbioportal.org/api'
studyId = 'luad_tcga_pan_can_atlas_2018'
geneSymbol = 'PTDSS1'

# Search cancer database for eligble studies 
def getStudies(): 
    response = requests.get(f'{BASE}/studies', timeout=30)
    response.raise_for_status()
    studies = response.json()

    for study in studies: 
        if "lung adenocarcinoma" in study["name"].lower():
            print(study["studyId"], study["name"])

## Pull the moleuculare profile of the study 
def getStudyMolecularProfile():
    response = requests.get(f'{BASE}/studies/{studyId}/molecular-profiles', timeout=30)
    response.raise_for_status()
    profiles = response.json()
    
    for profile in profiles: 
        if(profile['molecularAlterationType'] == 'COPY_NUMBER_ALTERATION' and profile['datatype'].lower() == 'discrete'):
              return(profile["molecularProfileId"])

## identifies if the study gas the clinical data I want to compare 
def getClinicalAttributes(): 
    response = requests.get(f'{BASE}/studies/{studyId}/clinical-attributes', timeout=30)
    response.raise_for_status()
    clinicalData = response.json()

    for data in clinicalData: 
        if data["clinicalAttributeId"] == "OS_MONTHS" or data["clinicalAttributeId"] == "OS_STATUS":
            print(data["clinicalAttributeId"])

## gets gene ID 
def getGeneInfo(): 
    response = requests.get(f'{BASE}/genes/{geneSymbol}', timeout=30)
    response.raise_for_status()
    geneInfo = response.json()

    return geneInfo["entrezGeneId"]

##gets sample
def getSampleList(): 
    response = requests.get(f'{BASE}/studies/{studyId}/sample-lists', timeout=30)
    response.raise_for_status()
    sampleList = response.json()

    for x in sampleList:
        if x['sampleListId'] == 'luad_tcga_pan_can_atlas_2018_cna':
            return x["sampleListId"]

##gets Id of samples from list 
def getSampleId():
    sampleListId = getSampleList()
    response = requests.get(f'{BASE}/sample-lists/{sampleListId}/sample-ids', timeout=30)
    response.raise_for_status()
    sampleIds = response.json()

    return(sampleIds[0:])
## Get the copy number of the sample 
def getCopyNumber():
    sampleIds = getSampleId()
    gene_id = getGeneInfo()
    molecularProfileId = getStudyMolecularProfile()

    body = {
        "sampleIds": sampleIds,
        "entrezGeneIds": [gene_id]
    }

    response = requests.post(f'{BASE}/molecular-profiles/{molecularProfileId}/discrete-copy-number/fetch', json=body,timeout=30, params={"discreteCopyNumberEventType": "ALL"})
    response.raise_for_status()
    copyNumber = response.json()

 
    return copyNumber

## Find the clinical data, in this case survival of patient 
def getClinicalData(patientId): 
    survival = {}
    
    response = requests.get(f'{BASE}/studies/{studyId}/patients/{patientId}/clinical-data')
    response.raise_for_status()
    patientData = response.json()

    for x in patientData:
        if x["clinicalAttributeId"] == 'OS_MONTHS' or x["clinicalAttributeId"] == "OS_STATUS":
            survival[x["clinicalAttributeId"]] = x["value"]

    return survival

def copyNumberAnalysis():
    rows = []
    copy_number_data = getCopyNumber()
    for record in copy_number_data: 
        survival = getClinicalData(record["patientId"])
        print(record["patientId"], record["alteration"], survival)

# print(getStudyMolecularProfile())
# getClinicalAttributes()
# print(getGeneInfo())
# print(getSampleList())
# print(getSampleId())
# getCopyNumber()
copyNumberAnalysis()
