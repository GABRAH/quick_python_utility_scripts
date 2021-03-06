import os
import shutil
import sys
import csv
import requests
from privateer import privateer_core as pvtcore
from privateer import privateer_modelling as pvtmodelling
from Bio import pairwise2
from Bio.Seq import Seq


def CreateFolder(path):
    if not os.path.exists(path):
        os.makedirs(path)


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


def get_sequences_in_receiving_model(receiverpath):
    builder_sequence_only = pvtmodelling.Builder(receiverpath, True)
    receiver_sequence = builder_sequence_only.get_receiving_model_sequence_info()

    return receiver_sequence


def import_uniprotID_list(path):
    print(path)

    fields = []
    output = []
    with open(path, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)

        for row in csvreader:
            # print(row)
            zip_iterator = zip(fields, row)
            zipped_dict = dict(zip_iterator)
            # print(zipped_dict)
            if not any(
                d.get("UniProtID", None) == zipped_dict["Uniprot ID"] for d in output
            ):
                entryToAppend = {
                    "UniProtID": zipped_dict["Uniprot ID"],
                    "glycosites": [int(zipped_dict["N-glycosylation site"])],
                }
                output.append(entryToAppend)

            if any(
                d.get("UniProtID", None) == zipped_dict["Uniprot ID"] for d in output
            ):
                index = next(
                    (
                        i
                        for i, item in enumerate(output)
                        if item["UniProtID"] == zipped_dict["Uniprot ID"]
                    ),
                    None,
                )
                if (
                    int(zipped_dict["N-glycosylation site"])
                    not in output[index]["glycosites"]
                ):
                    output[index]["glycosites"].append(
                        int(zipped_dict["N-glycosylation site"])
                    )

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
    PDBentries = [element for element in dbReferences if element["type"] == "PDB"]
    PDBentries_NMR = [
        element for element in PDBentries if element["properties"]["method"] == "NMR"
    ]
    PDBentries_XrayCryoEM = [
        element
        for element in PDBentries
        if element["properties"]["method"] == "X-ray"
        or element["properties"]["method"] == "EM"
    ]
    for index, item in enumerate(PDBentries_XrayCryoEM):
        try:
            oldString = item["properties"]["resolution"]
        except KeyError:
            PDBentries_XrayCryoEM[index]["properties"]["resolution"] = 42.0
            continue
        newString = oldString.replace(" A", "")
        convertedFloat = float(newString)
        PDBentries_XrayCryoEM[index]["properties"]["resolution"] = convertedFloat
    PDBentries_XrayCryoEM.sort(
        key=lambda e: e["properties"]["resolution"], reverse=False
    )
    for PDBentry in PDBentries_XrayCryoEM:
        chains = []
        modelID = PDBentry["id"]
        modelMethod = PDBentry["properties"]["method"]
        modelResolution = PDBentry["properties"]["resolution"]
        chainInfo = PDBentry["properties"]["chains"]
        # print(chainInfo)
        modelChains = chainInfo.split(",")
        for chain in modelChains:
            modelChainIDs = chain.split("=", 1)[0]
            modelChainRange = chain.split("=", 1)[1]
            modelChainStart = int(modelChainRange.split("-")[0])
            modelChainEnd = int(modelChainRange.split("-")[1])
            chains.append(
                {
                    "chainIDs": modelChainIDs,
                    "chain_start": modelChainStart,
                    "chain_end": modelChainEnd,
                }
            )
        # print(chains)
        entryToAppend = {
            "id": modelID,
            "method": modelMethod,
            "resolution": modelResolution,
            "chains": chains,
        }
        sanitised_PDBentries.append(entryToAppend)

    for PDBentry in PDBentries_NMR:
        chains = []
        modelID = PDBentry["id"]
        modelMethod = PDBentry["properties"]["method"]
        modelResolution = 0
        chainInfo = PDBentry["properties"]["chains"]
        # print(chainInfo)
        modelChains = chainInfo.split(",")
        for chain in modelChains:
            modelChainIDs = chain.split("=", 1)[0]
            modelChainRange = chain.split("=", 1)[1]
            modelChainStart = int(modelChainRange.split("-")[0])
            modelChainEnd = int(modelChainRange.split("-")[1])
            chains.append(
                {
                    "chainIDs": modelChainIDs,
                    "chain_start": modelChainStart,
                    "chain_end": modelChainEnd,
                }
            )
        # print(chains)
        entryToAppend = {
            "id": modelID,
            "method": modelMethod,
            "resolution": modelResolution,
            "chains": chains,
        }
        sanitised_PDBentries.append(entryToAppend)

    output = {
        "sequenceLength": outputSequenceLength,
        "sequence": outputSequence,
        "PDBs": sanitised_PDBentries,
    }

    return output


def query_PDBeKB(uniprotID, glycosites):
    sanitised_PDBentries = []
    PDBeKB_URL = f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/unipdb/{uniprotID}"
    PDBeKB_Response = requests.get(PDBeKB_URL, headers={"Accept": "application/json"})
    if not PDBeKB_Response.ok:
        return {"status": "failed"}
    PDBeKB_ResponseJSON = PDBeKB_Response.json()
    sequence = PDBeKB_ResponseJSON[uniprotID]["sequence"]
    experimentalData = PDBeKB_ResponseJSON[uniprotID]["data"]
    experimentalData.sort(
        key=lambda e: e["additionalData"]["resolution"], reverse=False
    )
    for glycosite in glycosites:
        PDBentries = []
        for entry in experimentalData:
            accessionID = entry["accession"]
            bestChainID = entry["bestChainId"]
            additionalData = entry["additionalData"]
            resolution = additionalData["resolution"]
            experiment = additionalData["experiment"]
            residues = entry["residues"]
            for residueRange in residues:
                start = residueRange["startIndex"]
                end = residueRange["endIndex"]
                glycositePresent = True if start <= glycosite <= end else False
                if glycositePresent == True:
                    break
            entryToAppend = {
                "PDBid": accessionID,
                "method": experiment,
                "resolution": resolution,
                "bestChain": bestChainID,
                "glycosite": glycosite,
                "glycositePresent": glycositePresent,
                "sequence": sequence,
            }
            PDBentries.append(entryToAppend)
        sanitised_PDBentries.append({"status": "success", "PDB": PDBentries})

    # print(sanitised_PDBentries)
    return sanitised_PDBentries


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
        entryToAppend = {
            "glycosite_present": glycosite_present,
            "chainID": chainID,
            "chainSequence": chainSequence,
        }
        output.append(entryToAppend)
    return output


def check_glycan_status_in_PDB_model(PDBentry, glycosite, outputDirectory):
    PDBID = PDBentry["PDBid"]
    tmpOutputDirectory = os.path.join(outputDirectory, "tmp")
    PDBIDFullPath = os.path.join(tmpOutputDirectory, f"{PDBID}.pdb")

    if not os.path.exists(PDBIDFullPath):
        if os.path.exists(tmpOutputDirectory):
            shutil.rmtree(tmpOutputDirectory)
        CreateFolder(tmpOutputDirectory)
        download_RCSBPDB_file(PDBID, tmpOutputDirectory)

    if os.path.exists(PDBIDFullPath):
        model_sequences = get_sequences_in_receiving_model(PDBIDFullPath)
        best_model_sequence = (
            PDBentry["bestChain"][0]
            if len(PDBentry["bestChain"]) > 1
            else PDBentry["bestChain"]
        )
        relevant_model_sequence = next(
            (
                item
                for item in model_sequences
                if item["ChainID"] == best_model_sequence
            ),
            None,
        )
        if relevant_model_sequence is not None:
            chain_residues = relevant_model_sequence["Residues"]
            uniprot_sequence = Seq(PDBentry["sequence"])
            pdb_sequence = Seq(relevant_model_sequence["Sequence"])
            alignments = pairwise2.align.localxx(pdb_sequence, uniprot_sequence)
            actual_alignment = alignments[0]
            zero_based_glycosite_index = glycosite - 1
            if (
                actual_alignment[0][zero_based_glycosite_index]
                == actual_alignment[1][zero_based_glycosite_index]
            ):
                recovered_sequence = []
                for i in range(glycosite):
                    char = actual_alignment[0][i]
                    if char.isalpha():
                        recovered_sequence.append(char)

                MiniMolIndex = len(recovered_sequence) - 1
                relevant_model_residue = next(
                    (item for item in chain_residues if item["index"] == MiniMolIndex),
                    None,
                )
                glycosylation = pvtcore.GlycosylationComposition_memsafe(PDBIDFullPath)
                glycanSummary = glycosylation.get_summary_of_detected_glycans()
                if len(glycanSummary):
                    relevant_glycan = next(
                        (
                            item
                            for item in glycanSummary
                            if item["RootInfo"]["ProteinChainID"] == best_model_sequence
                            and item["RootInfo"]["ProteinResidueSeqnum"]
                            == relevant_model_residue["residueSeqnum"]
                            and item["RootInfo"]["ProteinResidueType"]
                            == relevant_model_residue["residueType"]
                        ),
                        None,
                    )
                    if relevant_glycan is not None:
                        glycanWURCS = relevant_glycan["WURCS"]
                        glycanLength = glycanWURCS[12]
                        rootInfo = f'{relevant_glycan["RootInfo"]["ProteinResidueType"]}-{relevant_glycan["RootInfo"]["ProteinResidueSeqnum"]}'
                        output = {
                            "glycanLength": glycanLength,
                            "WURCS": glycanWURCS,
                            "rootInfo": rootInfo,
                        }
                        return output
                    else:
                        return None
                else:
                    return None
            else:
                return None
        else:
            return None
    else:
        return None


def get_pdbs_associated_with_uniprotid(
    uniprotList, outputDirectory, PDBs_containing_glycans
):
    output = []
    redirectedUniProtIDs = []
    master_csv_filename = "processed_uniprotIDs.csv"

    master_csv_fieldnames = [
        "original_UniProtID",
        "redirected_UniProtID",
        "glycosite_location",
        "PDB",
        "PDB_contains_glycans",
        "PDB_method",
        "resolution",
        "bestChain",
        "glycosite_present",
        "AlphaFold_confidence",
        "glycanLength",
        "glycosylated_PDB_amino_acid",
        "WURCS",
        "UniProt_sequence",
    ]
    PrepareFolderAndFile(outputDirectory, master_csv_filename)
    with open(
        os.path.join(outputDirectory, master_csv_filename), "a", newline=""
    ) as primaryresultsFile:
        primaryresultswriter = InitializeResultsFile(
            primaryresultsFile, master_csv_fieldnames
        )
        for index, item in enumerate(uniprotList):
            currentUniProtID = item["UniProtID"]
            currentGlycosites = item["glycosites"]
            print(f"Processing {index}/{len(uniprotList)} - {currentUniProtID}")
            try:
                uniprotQuery = query_uniprotID(currentUniProtID)
                redirected_uniprotID = None
                confidence = get_confidence_scores_from_alphafoldDB_model(
                    currentUniProtID, currentGlycosites
                )
            except requests.HTTPError as exception:
                redirected_uniprotID = get_redirect_link(currentUniProtID)
                uniprotQuery = query_uniprotID(redirected_uniprotID)
                confidence = get_confidence_scores_from_alphafoldDB_model(
                    redirected_uniprotID, currentGlycosites
                )
                redirectedUniProtIDs.append(
                    {
                        "old_UniProtID": currentUniProtID,
                        "old_UniProtID": redirected_uniprotID,
                    }
                )

            if redirected_uniprotID is None:
                PDBeKB_query = query_PDBeKB(currentUniProtID, currentGlycosites)
            else:
                PDBeKB_query = query_PDBeKB(redirected_uniprotID, currentGlycosites)

            for glycositeIndex, glycosite in enumerate(currentGlycosites):
                glycositeConfidenceDict = next(
                    (e for e in confidence if e["glycosite"] == glycosite), None
                )
                if glycositeConfidenceDict is not None:
                    glycositeConfidence = glycositeConfidenceDict["confidence"]
                else:
                    glycositeConfidence = "-"

                try:
                    currentGlycositePDBeKBQuery = PDBeKB_query[glycositeIndex]
                    if currentGlycositePDBeKBQuery["status"] == "success":
                        if len(currentGlycositePDBeKBQuery["PDB"]):
                            for PDBentry in currentGlycositePDBeKBQuery["PDB"]:
                                glycanStatus = check_glycan_status_in_PDB_model(
                                    PDBentry, glycosite, outputDirectory
                                )
                                if glycanStatus is not None:
                                    currentRow = {
                                        "original_UniProtID": currentUniProtID,
                                        "redirected_UniProtID": "-"
                                        if redirected_uniprotID is None
                                        else redirected_uniprotID,
                                        "glycosite_location": glycosite,
                                        "PDB": PDBentry["PDBid"],
                                        "PDB_contains_glycans": "No"
                                        if not any(
                                            dict.get("PDBID", "No") == PDBentry["PDBid"]
                                            for dict in PDBs_containing_glycans
                                        )
                                        else "Yes",
                                        "PDB_method": PDBentry["method"],
                                        "resolution": PDBentry["resolution"],
                                        "bestChain": PDBentry["bestChain"],
                                        "glycosite_present": "Yes"
                                        if PDBentry["glycositePresent"] == True
                                        else "No",
                                        "AlphaFold_confidence": glycositeConfidence,
                                        "glycanLength": glycanStatus["glycanLength"],
                                        "glycosylated_PDB_amino_acid": glycanStatus[
                                            "rootInfo"
                                        ],
                                        "WURCS": glycanStatus["WURCS"],
                                        "UniProt_sequence": PDBentry["sequence"],
                                    }
                                    primaryresultswriter.writerow(currentRow)
                                    output.append(currentRow)
                                else:
                                    currentRow = {
                                        "original_UniProtID": currentUniProtID,
                                        "redirected_UniProtID": "-"
                                        if redirected_uniprotID is None
                                        else redirected_uniprotID,
                                        "glycosite_location": glycosite,
                                        "PDB": PDBentry["PDBid"],
                                        "PDB_contains_glycans": "No"
                                        if not any(
                                            dict.get("PDBID", "No") == PDBentry["PDBid"]
                                            for dict in PDBs_containing_glycans
                                        )
                                        else "Yes",
                                        "PDB_method": PDBentry["method"],
                                        "resolution": PDBentry["resolution"],
                                        "bestChain": PDBentry["bestChain"],
                                        "glycosite_present": "Yes"
                                        if PDBentry["glycositePresent"] == True
                                        else "No",
                                        "AlphaFold_confidence": glycositeConfidence,
                                        "glycanLength": "-",
                                        "glycosylated_PDB_amino_acid": "-",
                                        "WURCS": "-",
                                        "UniProt_sequence": PDBentry["sequence"],
                                    }
                                    primaryresultswriter.writerow(currentRow)
                                    output.append(currentRow)
                        else:
                            currentRow = {
                                "original_UniProtID": currentUniProtID,
                                "redirected_UniProtID": "-"
                                if redirected_uniprotID is None
                                else redirected_uniprotID,
                                "glycosite_location": glycosite,
                                "PDB": "-",
                                "PDB_contains_glycans": "-",
                                "PDB_method": "-",
                                "resolution": "-",
                                "bestChain": "-",
                                "glycosite_present": "-",
                                "AlphaFold_confidence": glycositeConfidence,
                                "glycanLength": "-",
                                "glycosylated_PDB_amino_acid": "-",
                                "WURCS": "-",
                                "UniProt_sequence": "-",
                            }
                            primaryresultswriter.writerow(currentRow)
                            output.append(currentRow)
                    else:
                        currentRow = {
                            "original_UniProtID": currentUniProtID,
                            "redirected_UniProtID": "-"
                            if redirected_uniprotID is None
                            else redirected_uniprotID,
                            "glycosite_location": glycosite,
                            "PDB": "-",
                            "PDB_contains_glycans": "-",
                            "PDB_method": "-",
                            "resolution": "-",
                            "bestChain": "-",
                            "glycosite_present": "-",
                            "AlphaFold_confidence": glycositeConfidence,
                            "glycanLength": "-",
                            "glycosylated_PDB_amino_acid": "-",
                            "WURCS": "-",
                            "UniProt_sequence": "-",
                        }
                        primaryresultswriter.writerow(currentRow)
                        output.append(currentRow)
                except (IndexError, KeyError) as e:
                    currentRow = {
                        "original_UniProtID": currentUniProtID,
                        "redirected_UniProtID": "-"
                        if redirected_uniprotID is None
                        else redirected_uniprotID,
                        "glycosite_location": glycosite,
                        "PDB": "-",
                        "PDB_contains_glycans": "-",
                        "PDB_method": "-",
                        "resolution": "-",
                        "bestChain": "-",
                        "glycosite_present": "-",
                        "AlphaFold_confidence": glycositeConfidence,
                        "glycanLength": "-",
                        "glycosylated_PDB_amino_acid": "-",
                        "WURCS": "-",
                        "UniProt_sequence": "-",
                    }
                    primaryresultswriter.writerow(currentRow)
                    output.append(currentRow)

    return output


def get_confidence_scores_from_alphafoldDB_model(uniprotID, glycosites):
    output = []
    savedLines = []
    try:
        requestURL = f"https://alphafold.ebi.ac.uk/files/AF-{uniprotID}-F1-model_v1.pdb"
        query = requests.get(requestURL, allow_redirects=True)

        outputLines = []
        downloadedLines = query.iter_lines()
        for line in downloadedLines:
            savedLines.append(line)
    except Exception as e:
        print(f"Error: {e}")

    for glycosite in glycosites:
        for line in savedLines:
            decodedLine = line.decode("utf-8")
            if decodedLine[:4] == "ATOM":
                residueNumber = int(decodedLine[22:26])
                if residueNumber == glycosite:
                    confidenceNumber = float(decodedLine[61:66])
                    output.append(
                        {"glycosite": glycosite, "confidence": confidenceNumber}
                    )
                    break

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
    failedDownload = "<p>The requested URL was not found on this server.</p>"
    failure = False
    for line in downloadedLines:
        decodedLine = line.decode("utf-8")
        if failedDownload in decodedLine:
            failure = True
            break
        else:
            outputLines.append(decodedLine)

    if failure == False:
        with open(outputFilePath, "w") as file:
            file.writelines("%s\n" % l for l in outputLines)
        print(
            f"\tSuccessfully downloaded model from RCSB PDB with PDB ID: '{PDBID}' to {outputFilePath}"
        )
    else:
        print(
            f"\tFailed to download model from RCSB PDB with PDB ID: '{PDBID}' to {outputFilePath}"
        )
    return True


def download_models_for_glycosite_view(csvoutput, outputFolderLocation):
    glycositeViewFolder = os.path.join(outputFolderLocation, "glycositeViewAndModels")
    for index, row in enumerate(csvoutput):
        print(
            f'Looking to download for {index}/{len(csvoutput)} - {row["original_UniProtID"]}'
        )
        if row["glycosite_present"] == "Yes":
            print(
                f'\tActually downloading PDB for {index}/{len(csvoutput)} - {row["original_UniProtID"]}'
            )
            uniProtID = (
                row["original_UniProtID"]
                if row["redirected_UniProtID"] == "-"
                else f'{row["original_UniProtID"]}__{row["redirected_UniProtID"]}'
            )
            glycosite = row["glycosite_location"]
            PDBID = row["PDB"]
            newFolder = os.path.join(glycositeViewFolder, f"{uniProtID}/{glycosite}")
            CreateFolder(newFolder)
            download_RCSBPDB_file(PDBID, newFolder)
            download_and_prepare_alphafoldDB_model(row["original_UniProtID"], newFolder)


def import_generic_csv_file(path):
    print(path)
    output = []
    with open(os.path.join(path)) as previous_output:
        imported_dict = csv.DictReader(previous_output)
        for row in imported_dict:
            output.append(row)

    return output


def get_list_of_PDB_containing_glycans(uniprotOutput, privateerOutput):
    output = []

    for row_uniprot in uniprotOutput:
        uniprotID = row_uniprot["original_UniProtID"]
        glycosite_present = row_uniprot["glycosite_present"]
        PDBID = row_uniprot["PDB"]
        if glycosite_present == "Yes":
            matching_entry = next(
                (item for item in privateerOutput if item["pdbID"] == PDBID),
                None,
            )
            if matching_entry is not None:
                already_appended = next(
                    (item for item in output if item["PDBID"] == PDBID),
                    None,
                )
                if already_appended is None:
                    output.append(
                        {"uniprotID": uniprotID, "PDBID": matching_entry["pdbID"]}
                    )

    return output


if __name__ == "__main__":
    fileLocation = os.path.abspath(__file__)
    nDirectoriesToGoBack = 1
    projectDir = os.path.normpath(
        os.path.join(*([fileLocation] + [".."] * nDirectoriesToGoBack))
    )
    inputPathFileLocationAlpha = os.path.join(projectDir, "input/uniprotIDs.csv")
    inputPathFileLocationBravo = os.path.join(
        projectDir, "input/previously_processed_uniprotIDs.csv"
    )
    inputPathFileLocationCharlie = os.path.join(
        projectDir, "input/processed_pdbmirror.csv"
    )
    outputFolderLocation = os.path.join(projectDir, "output")
    PrepareFolder(outputFolderLocation)
    previousUniProtOutput = import_generic_csv_file(inputPathFileLocationBravo)
    privateerCurationOutput = import_generic_csv_file(inputPathFileLocationCharlie)
    PDBs_containing_glycans = get_list_of_PDB_containing_glycans(
        previousUniProtOutput, privateerCurationOutput
    )
    uniprotList = import_uniprotID_list(inputPathFileLocationAlpha)
    outputted = get_pdbs_associated_with_uniprotid(
        uniprotList, outputFolderLocation, PDBs_containing_glycans
    )
    download_models_for_glycosite_view(outputted, outputFolderLocation)

    print("Script has finished running!")
