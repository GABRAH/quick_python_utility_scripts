import os
import shutil
import sys
import csv
import requests

def PrepareFolder(path):
    if os.path.exists(path):
        shutil.rmtree(path)
        os.makedirs(path)
    if not os.path.exists(path):
        os.makedirs(path)

def PrepareFolderAndFile(path, outputFileName):
    if not os.path.exists(path):
        os.makedirs(path)
    if os.path.exists(os.path.join(path, outputFileName)):
        os.remove(os.path.join(path, outputFileName))

def InitializeResultsFile(inputFile, fieldnames):
    writer = csv.DictWriter(inputFile, fieldnames=fieldnames)
    currentRow = dict.fromkeys(fieldnames, None)
    for key, value in currentRow.items():
        if value is None:
            currentRow[key] = key
    writer.writerow(currentRow)
    return writer


def import_uniprotID_list(path):
    print(path)
    
    fields = []
    output = []
    with open(path, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)
    
        for row in csvreader:
            # print(row)
            zip_iterator = zip(fields, row)
            zipped_dict = dict(zip_iterator)
            # print(zipped_dict)
            if not any(d.get('UniProtID', None) == zipped_dict["Uniprot ID"] for d in output):
                entryToAppend = {"UniProtID": zipped_dict["Uniprot ID"], "glycosites": [int(zipped_dict["N-glycosylation site"])]}
                output.append(entryToAppend)
            
            if any(d.get('UniProtID', None) == zipped_dict["Uniprot ID"] for d in output):
                index = next((i for i, item in enumerate(output) if item["UniProtID"] == zipped_dict["Uniprot ID"]), None)
                if int(zipped_dict["N-glycosylation site"]) not in output[index]["glycosites"]:
                    output[index]["glycosites"].append(int(zipped_dict["N-glycosylation site"]))
            
    return output


def query_uniprotID(uniprotID):
    sanitised_PDBentries = []
    uniprotRequestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprotID}"
    uniprotResponse = requests.get(
        uniprotRequestURL, headers={"Accept": "application/json"}
    )
    if not uniprotResponse.ok:
        uniprotResponse.raise_for_status()
        sys.exit()
    uniprotResponseJSON = uniprotResponse.json()
    uniprotSequence = uniprotResponseJSON["sequence"]
    outputSequence = uniprotSequence["sequence"]
    outputSequenceLength = uniprotSequence["length"]
    dbReferences = uniprotResponseJSON["dbReferences"]
    PDBentries = [element for element in dbReferences if element['type'] == "PDB"]
    PDBentries_NMR = [element for element in PDBentries if element['properties']['method'] == "NMR"]
    PDBentries_XrayCryoEM = [element for element in PDBentries if element['properties']['method'] == "X-ray" or element['properties']['method'] == "EM"]
    for index, item in enumerate(PDBentries_XrayCryoEM):
        try:
            oldString = item['properties']['resolution']
        except KeyError:
            PDBentries_XrayCryoEM[index]['properties']['resolution'] = 42.0
            continue
        newString = oldString.replace(' A', '')
        convertedFloat = float(newString)
        PDBentries_XrayCryoEM[index]['properties']['resolution'] = convertedFloat
    PDBentries_XrayCryoEM.sort(key=lambda e: e['properties']['resolution'], reverse=False)
    for PDBentry in PDBentries_XrayCryoEM:
        chains = []
        modelID = PDBentry['id']
        modelMethod = PDBentry['properties']['method']
        modelResolution = PDBentry['properties']['resolution']
        chainInfo = PDBentry['properties']['chains']
        # print(chainInfo)
        modelChains = chainInfo.split(",")
        for chain in modelChains:
            modelChainIDs = chain.split("=", 1)[0]
            modelChainRange = chain.split("=", 1)[1]
            modelChainStart = int(modelChainRange.split("-")[0])
            modelChainEnd = int(modelChainRange.split("-")[1])
            chains.append({"chainIDs": modelChainIDs, "chain_start": modelChainStart, "chain_end": modelChainEnd})
        # print(chains)
        entryToAppend = {"id": modelID, "method": modelMethod, "resolution": modelResolution, "chains": chains}
        sanitised_PDBentries.append(entryToAppend)
    
    for PDBentry in PDBentries_NMR:
        chains = []
        modelID = PDBentry['id']
        modelMethod = PDBentry['properties']['method']
        modelResolution = 0
        chainInfo = PDBentry['properties']['chains']
        # print(chainInfo)
        modelChains = chainInfo.split(",")
        for chain in modelChains:
            modelChainIDs = chain.split("=", 1)[0]
            modelChainRange = chain.split("=", 1)[1]
            modelChainStart = int(modelChainRange.split("-")[0])
            modelChainEnd = int(modelChainRange.split("-")[1])
            chains.append({"chainIDs": modelChainIDs, "chain_start": modelChainStart, "chain_end": modelChainEnd})
        # print(chains)
        entryToAppend = {"id": modelID, "method": modelMethod, "resolution": modelResolution, "chains": chains}
        sanitised_PDBentries.append(entryToAppend)

    output = {
        "sequenceLength": outputSequenceLength,
        "sequence": outputSequence,
        "PDBs": sanitised_PDBentries
    }
        
    return output

def query_PDBeKB(uniprotID):
    pass


def get_redirect_link(uniprotID):
    r = requests.get(f"https://www.uniprot.org/uniprot/{uniprotID}")
    redirectedURL = r.url
    splitURL = redirectedURL.split("/")
    redirectedUniProtID = splitURL[-1]
    return redirectedUniProtID

def check_if_glycosite_present_in_PDB(glycosite, chains, uniprot_sequence):
    output = []
    for chain in chains:
        chainSequence = ""
        chain_start = chain["chain_start"]
        chain_end = chain["chain_end"]
        chainID = chain["chainIDs"]
        glycosite_present = "Yes" if chain_start <= glycosite <= chain_end else "No"
        zero_index_based_chain_start = 0 if chain_start < 1 else chain_start - 1
        for index in range(zero_index_based_chain_start, chain_end):
            chainSequence += uniprot_sequence[index]
        entryToAppend = {"glycosite_present": glycosite_present, "chainID": chainID, "chainSequence": chainSequence}
        output.append(entryToAppend)
    return output

            
        

def get_pdbs_associated_with_uniprotid(uniprotList, outputDirectory):
    output = []
    redirectedUniProtIDs = []
    master_csv_filename = "processed_uniprotIDs.csv"
    
    master_csv_fieldnames = [
        "original_UniProtID",
        "redirected_UniProtID",
        "glycosite_location",
        "PDB",
        "PDB_method",
        "resolution",
        "glycosite_present_in_PDB",
        "glycosite_present_in_PDB_ChainID",
        "UniProt_sequence",
        "PDB_ChainID_1",
        "PDB_sequence_1", 
        "PDB_ChainID_2",
        "PDB_sequence_2",
        "PDB_ChainID_3",
        "PDB_sequence_3"
    ]
    PrepareFolderAndFile(outputDirectory, master_csv_filename)
    with open(
        os.path.join(outputDirectory, master_csv_filename), "a", newline=""
    ) as primaryresultsFile:
        primaryresultswriter = InitializeResultsFile(primaryresultsFile, master_csv_fieldnames)
        for index, item in enumerate(uniprotList):
            currentUniProtID = item["UniProtID"]
            currentGlycosites = item["glycosites"]
            print(f'Processing {index}/{len(uniprotList)} - {currentUniProtID}')
            try:
                uniprotQuery = query_uniprotID(currentUniProtID)
                redirected_uniprotID = None
            except requests.HTTPError as exception:
                redirected_uniprotID = get_redirect_link(currentUniProtID)
                uniprotQuery = query_uniprotID(redirected_uniprotID)
                redirectedUniProtIDs.append({"old_UniProtID": currentUniProtID, "old_UniProtID": redirected_uniprotID})
        
            associatedPDBs = uniprotQuery["PDBs"]
            for glycosite in currentGlycosites:
                numberOfAssociatedPDBs = len(associatedPDBs)
                if numberOfAssociatedPDBs == 0:
                    currentRow = {
                        "original_UniProtID": currentUniProtID,
                        "redirected_UniProtID": "-" if redirected_uniprotID is None else redirected_uniprotID,
                        "glycosite_location": glycosite,
                        "PDB": "-",
                        "PDB_method": "-",
                        "resolution": "-",
                        "glycosite_present_in_PDB": "-",
                        "glycosite_present_in_PDB_ChainID": "-",
                        "UniProt_sequence": uniprotQuery["sequence"],
                        "PDB_ChainID_1": "-",
                        "PDB_sequence_1": "-", 
                        "PDB_ChainID_2": "-",
                        "PDB_sequence_2": "-",
                        "PDB_ChainID_3": "-",
                        "PDB_sequence_3": "-"
                    }
                    primaryresultswriter.writerow(currentRow)
                    output.append(currentRow)
                else:
                    for associatedPDB in associatedPDBs:
                        numberOfChains = len(associatedPDB["chains"])
                        if numberOfChains == 1:
                            chains = associatedPDB["chains"]
                            presence_list = check_if_glycosite_present_in_PDB(glycosite, chains, uniprotQuery["sequence"])
                            currentRow = {
                                "original_UniProtID": currentUniProtID,
                                "redirected_UniProtID": "-" if redirected_uniprotID is None else redirected_uniprotID,
                                "glycosite_location": glycosite,
                                "PDB": associatedPDB["id"],
                                "PDB_method": associatedPDB["method"],
                                "resolution": associatedPDB["resolution"],
                                "glycosite_present_in_PDB": presence_list[0]["glycosite_present"],
                                "glycosite_present_in_PDB_ChainID": presence_list[0]["chainID"],
                                "UniProt_sequence": uniprotQuery["sequence"],
                                "PDB_ChainID_1": presence_list[0]["chainID"],
                                "PDB_sequence_1": presence_list[0]["chainSequence"], 
                                "PDB_ChainID_2": "-",
                                "PDB_sequence_2": "-",
                                "PDB_ChainID_3": "-",
                                "PDB_sequence_3": "-"
                            }
                            primaryresultswriter.writerow(currentRow)
                            output.append(currentRow)
                        elif numberOfChains == 2:
                            chains = associatedPDB["chains"]
                            presence_list = check_if_glycosite_present_in_PDB(glycosite, chains, uniprotQuery["sequence"])
                            glycosite_present_in_PDB = "No"
                            glycosite_present_in_chains = "-"
                            for index, item in enumerate(presence_list):
                                if item["glycosite_present"] == "Yes":
                                    glycosite_present_in_PDB = "Yes"
                                    glycosite_present_in_chains = item["chainID"]
                                    glycosite_present_in_chains += "/" if index < len(presence_list) else ""
                            
                            currentRow = {
                                "original_UniProtID": currentUniProtID,
                                "redirected_UniProtID": "-" if redirected_uniprotID is None else redirected_uniprotID,
                                "glycosite_location": glycosite,
                                "PDB": associatedPDB["id"],
                                "PDB_method": associatedPDB["method"],
                                "resolution": associatedPDB["resolution"],
                                "glycosite_present_in_PDB": glycosite_present_in_PDB,
                                "glycosite_present_in_PDB_ChainID": glycosite_present_in_chains,
                                "UniProt_sequence": uniprotQuery["sequence"],
                                "PDB_ChainID_1": presence_list[0]["chainID"],
                                "PDB_sequence_1": presence_list[0]["chainSequence"], 
                                "PDB_ChainID_2": presence_list[1]["chainID"],
                                "PDB_sequence_2": presence_list[1]["chainSequence"],
                                "PDB_ChainID_3": "-",
                                "PDB_sequence_3": "-"
                            }
                            primaryresultswriter.writerow(currentRow)
                            output.append(currentRow)
                        elif numberOfChains == 3:
                            chains = associatedPDB["chains"]
                            presence_list = check_if_glycosite_present_in_PDB(glycosite, chains, uniprotQuery["sequence"])
                            glycosite_present_in_PDB = "No"
                            glycosite_present_in_chains = "-"
                            for index, item in enumerate(presence_list):
                                if item["glycosite_present"] == "Yes":
                                    glycosite_present_in_PDB = "Yes"
                                    glycosite_present_in_chains = item["chainID"]
                                    glycosite_present_in_chains += "/" if index < len(presence_list) else ""
                            
                            currentRow = {
                                "original_UniProtID": currentUniProtID,
                                "redirected_UniProtID": "-" if redirected_uniprotID is None else redirected_uniprotID,
                                "glycosite_location": glycosite,
                                "PDB": associatedPDB["id"],
                                "PDB_method": associatedPDB["method"],
                                "resolution": associatedPDB["resolution"],
                                "glycosite_present_in_PDB": glycosite_present_in_PDB,
                                "glycosite_present_in_PDB_ChainID": glycosite_present_in_chains,
                                "UniProt_sequence": uniprotQuery["sequence"],
                                "PDB_ChainID_1": presence_list[0]["chainID"],
                                "PDB_sequence_1": presence_list[0]["chainSequence"], 
                                "PDB_ChainID_2": presence_list[1]["chainID"],
                                "PDB_sequence_2": presence_list[1]["chainSequence"],
                                "PDB_ChainID_3": presence_list[2]["chainID"],
                                "PDB_sequence_3": presence_list[2]["chainSequence"]
                            }
                            primaryresultswriter.writerow(currentRow)
                            output.append(currentRow)
                        else:
                            raise Exception(f"There are more chains than anticipated, expecting 3, got {numberOfChains}")
    return output
            


        

def download_and_prepare_alphafoldDB_model(uniprotID, downloadLocation):
    outputFileName = "AlphaFold_" + uniprotID + ".pdb"
    outputFilePath = os.path.join(downloadLocation, outputFileName)
    try:
        requestURL = f"https://alphafold.ebi.ac.uk/files/AF-{uniprotID}-F1-model_v1.pdb"
        query = requests.get(requestURL, allow_redirects=True)

        outputLines = []
        downloadedLines = query.iter_lines()
        for line in downloadedLines:
            decodedLine = line.decode("utf-8")
            if decodedLine[:5] != "MODEL":
                outputLines.append(decodedLine)

        with open(outputFilePath, "w") as file:
            file.writelines("%s\n" % l for l in outputLines)

        print(
            f"\tSuccessfully downloaded model from AlphaFoldDB with UniProt ID: '{uniprotID}' to {outputFilePath}"
        )
    except Exception as e:
        print(f"Error: {e}")
        

def download_RCSBPDB_file(PDBID, downloadLocation):
    outputFileName = PDBID + ".pdb"
    outputFilePath = os.path.join(downloadLocation, outputFileName)
    requestURL = f"https://files.rcsb.org/download/{PDBID}.pdb"
    query = requests.get(requestURL, allow_redirects=True)

    outputLines = []
    downloadedLines = query.iter_lines()
    for line in downloadedLines:
        outputLines.append(line)
    
    with open(outputFilePath, "w") as file:
        file.writelines("%s\n" % l for l in outputLines)
    
    print(
        f"\tSuccessfully downloaded model from RCSB PDB with PDB ID: '{PDBID}' to {outputFilePath}"
    )


def download_models_for_glycosite_view(csvoutput, outputFolderLocation):
    glycositeViewFolder = os.path.join(outputFolderLocation, "glycositeViewAndModels")
    for index, row in enumerate(csvoutput):
        print(f'Looking to download for {index}/{len(csvoutput)} - {row["original_UniProtID"]}')
        if row["glycosite_present_in_PDB"] == "Yes":
            print(f'\tActually downloading PDB for {index}/{len(csvoutput)} - {row["original_UniProtID"]}')
            uniProtID = row["original_UniProtID"] if row["redirected_UniProtID"] == "-" else f'{row["original_UniProtID"]}__{row["redirected_UniProtID"]}'
            glycosite = row["glycosite_location"]
            PDBID = row["PDB"]
            newFolder = os.path.join(glycositeViewFolder, f'{uniProtID}/{glycosite}')
            PrepareFolder(newFolder)
            download_RCSBPDB_file(PDBID, newFolder)
            download_and_prepare_alphafoldDB_model(row["original_UniProtID"], newFolder)
            

            



if __name__ == "__main__":
    fileLocation = os.path.abspath(__file__)
    nDirectoriesToGoBack = 1
    projectDir = os.path.normpath(
        os.path.join(*([fileLocation] + [".."] * nDirectoriesToGoBack))
    )
    inputPathFileLocation = os.path.join(projectDir, "input/uniprotIDs.csv")
    outputFolderLocation = os.path.join(projectDir, "output")
    PrepareFolder(outputFolderLocation)
    uniprotList = import_uniprotID_list(inputPathFileLocation)
    # print(uniprotList)
    outputted = get_pdbs_associated_with_uniprotid(uniprotList, outputFolderLocation)
    download_models_for_glycosite_view(outputted, outputFolderLocation)
    print("Script has finished running!")
    