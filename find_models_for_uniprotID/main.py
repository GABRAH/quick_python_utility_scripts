import os
import re
import sys
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
    # need to convert here resolution field to float by removing amstrong symbol
    PDBentries_XrayCryoEM_sorted = PDBentries_XrayCryoEM.sort(key=lambda e: float(e['properties']['resolution']), reverse=True)
    for PDBentry in PDBentries_XrayCryoEM_sorted:
        print(PDBentry)
    
    print(outputSequence)
    print(outputSequenceLength)
    output = {
        "sequenceLength": outputSequenceLength,
        "sequence": outputSequence,
    }

    return "blabla"


if __name__ == "__main__":
    fileLocation = os.path.abspath(__file__)
    uniprotID = "P07602"
    query_uniprotID(uniprotID=uniprotID)