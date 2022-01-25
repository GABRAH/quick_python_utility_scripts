import os
import re
import sys
import csv
import argparse
import requests
import warnings
import json


def download_and_prepare_alphafoldDB_model(uniprotID, downloadLocation):
    outputFileName = uniprotID + ".pdb"
    outputFilePath = os.path.join(downloadLocation, outputFileName)
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
        f"Successfully downloaded model from AlphaFoldDB with UniProt ID: '{uniprotID}' to {outputFilePath}"
    )
    return outputFilePath


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
        oldString = item['properties']['resolution']
        newString = oldString.replace(' A', '')
        convertedFloat = float(newString)
        PDBentries_XrayCryoEM[index]['properties']['resolution'] = convertedFloat
    PDBentries_XrayCryoEM.sort(key=lambda e: e['properties']['resolution'], reverse=False)
    for PDBentry in PDBentries_XrayCryoEM:
        modelID = PDBentry['id']
        modelMethod = PDBentry['properties']['method']
        modelResolution = PDBentry['properties']['resolution']
        chainInfo = PDBentry['properties']['chains']
        modelChains = chainInfo.split("=", 1)[0]
        modelChainRange = chainInfo.split("=", 1)[1]
        modelChainStart = int(modelChainRange.split("-")[0])
        modelChainEnd = int(modelChainRange.split("-")[1])
        entryToAppend = {"id": modelID, "method": modelMethod, "resolution": modelResolution, "chains": modelChains, "chain_start": modelChainStart, "chain_end": modelChainEnd}
        sanitised_PDBentries.append(entryToAppend)
    
    for PDBentry in PDBentries_NMR:
        modelID = PDBentry['id']
        modelMethod = PDBentry['properties']['method']
        modelResolution = 0
        chainInfo = PDBentry['properties']['chains']
        modelChains = chainInfo.split("=", 1)[0]
        modelChainRange = chainInfo.split("=", 1)[1]
        modelChainStart = int(modelChainRange.split("-")[0])
        modelChainEnd = int(modelChainRange.split("-")[1])
        entryToAppend = {"id": modelID, "method": modelMethod, "resolution": modelResolution, "chains": modelChains, "chain_start": modelChainStart, "chain_end": modelChainEnd}
        sanitised_PDBentries.append(entryToAppend)

        output = {
            "sequenceLength": outputSequenceLength,
            "sequence": outputSequence,
            "PDBs": sanitised_PDBentries
        }
        
        return output

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
            if not any(d.get('UniprotID', None) == zipped_dict["Uniprot ID"] for d in output):
                entryToAppend = {"UniprotID": zipped_dict["Uniprot ID"], "glycosites": [int(zipped_dict["N-glycosylation site"])]}
                output.append(entryToAppend)
            
            if any(d.get('UniprotID', None) == zipped_dict["Uniprot ID"] for d in output):
                index = next((i for i, item in enumerate(output) if item["UniprotID"] == zipped_dict["Uniprot ID"]), None)
                if int(zipped_dict["N-glycosylation site"]) not in output[index]["glycosites"]:
                    output[index]["glycosites"].append(int(zipped_dict["N-glycosylation site"]))
            
    return output
            



if __name__ == "__main__":
    fileLocation = os.path.abspath(__file__)
    nDirectoriesToGoBack = 1
    projectDir = os.path.normpath(
        os.path.join(*([fileLocation] + [".."] * nDirectoriesToGoBack))
    )
    inputPathFileLocation = os.path.join(projectDir, "input/uniprotIDs.csv")
    outputFolderLocation = os.path.join(projectDir, "output")
    uniprotList = import_uniprotID_list(inputPathFileLocation)
    
    # get_pdbs_associated_with_uniprotid
    
    
    
    # uniprotID = "P07602"
    # query = query_uniprotID(uniprotID=uniprotID)
    