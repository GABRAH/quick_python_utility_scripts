import os, shutil
import csv
import json
from datetime import date

input_bins_relative_path = "input/bins"
output_directory = "output"


def CreateFolder(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def get_path_of_directory_file_is_located_in():
    directory = os.path.dirname(os.path.abspath(__file__))
    return directory


def parse_header_line(input_line):
    output = {"count_mean": None, "count_stdev": None}
    expected_num_of_values = 2
    current_string = input_line[0].strip()
    space_delimited_values = current_string.split(" ")
    if len(space_delimited_values
           ) == expected_num_of_values and "nan" not in space_delimited_values:
        count_mean = float(space_delimited_values[0])
        count_stdev = float(space_delimited_values[1])

        output["count_mean"] = count_mean
        output["count_stdev"] = count_stdev

        return output
    else:
        return output


def parse_current_nonheader_line_via_map(input_line):
    output = {
        "lower_phi": None,
        "higher_phi": None,
        "lower_psi": None,
        "higher_psi": None,
        "count": None
    }
    expected_num_of_values = 5
    current_string = input_line[0].strip()
    space_delimited_values = current_string.split(" ")
    if len(space_delimited_values
           ) == expected_num_of_values and "nan" not in space_delimited_values:
        lower_phi = float(space_delimited_values[0])
        higher_phi = float(space_delimited_values[1])
        lower_psi = float(space_delimited_values[2])
        higher_psi = float(space_delimited_values[3])
        count = float(space_delimited_values[4])

        output["lower_phi"] = lower_phi
        output["higher_phi"] = higher_phi
        output["lower_psi"] = lower_psi
        output["higher_psi"] = higher_psi
        output["count"] = count

        return output
    else:
        return output


def parse_nonheader_lines(input_lines):
    parsed_lines = list(map(parse_current_nonheader_line_via_map, input_lines))
    filtered_list = [
        current_dict for current_dict in parsed_lines
        if (current_dict["count"] > 0 and "nan" not in current_dict.values())
    ]
    return filtered_list


def process_current_dotdata_file(current_dotdata_file_path):
    # tmpoutput = {"path": current_dotdata_file_path}
    output = {"header": None, "bin_data": None}
    list_of_newline_delimited_rows = list(
        csv.reader(open(current_dotdata_file_path, 'r'), delimiter='\n'))

    header_line = parse_header_line(list_of_newline_delimited_rows[0])
    parsed_nonheader_list = parse_nonheader_lines(
        list_of_newline_delimited_rows[1:])

    output["header"] = header_line
    output["bin_data"] = parsed_nonheader_list
    return output


def parse_input_filename(filename):
    output = {
        "donor": None,
        "acceptor": None,
        "Linkage": {
            "donor_end": None,
            "acceptor_end": None
        }
    }

    sugar_info_only = filename.split(".")[0]
    strings_split_by_comma = sugar_info_only.split(",")
    acceptor_with_dash = strings_split_by_comma[1]
    donor_with_dash = strings_split_by_comma[0]

    output["donor"] = donor_with_dash.split("-")[0]
    output["acceptor"] = acceptor_with_dash.split("-")[1]
    output["Linkage"]["donor_end"] = acceptor_with_dash.split("-")[0]
    output["Linkage"]["acceptor_end"] = donor_with_dash.split("-")[1]

    print(f"filename: {filename} \t {output}")
    return output


def generate_output_for_unique_donor_sugar(sugar_info, current_linkage_data):
    output_dict = {
        "donor":
        sugar_info["donor"],
        "acceptor": [{
            "sugar":
            sugar_info["acceptor"],
            "Linkage": [{
                "donor_end": sugar_info["Linkage"]["donor_end"],
                "acceptor_end": sugar_info["Linkage"]["acceptor_end"],
                "Linkage_data": {
                    "summary": current_linkage_data["header"],
                    "bin_data": current_linkage_data["bin_data"]
                }
            }]
        }],
    }
    return output_dict


def append_output_for_unique_acceptor_sugar(iterator, acceptor_sugar, Linkage,
                                            current_linkage_data):
    iterator["acceptor"].append({
        "sugar":
        acceptor_sugar,
        "Linkage": [{
            "donor_end": Linkage["donor_end"],
            "acceptor_end": Linkage["acceptor_end"],
            "Linkage_data": {
                "summary": current_linkage_data["header"],
                "bin_data": current_linkage_data["bin_data"]
            }
        }]
    })


def append_output_for_unique_linkage(iterator, Linkage, current_linkage_data):
    iterator["Linkage"].append({
        "donor_end": Linkage["donor_end"],
        "acceptor_end": Linkage["acceptor_end"],
        "Linkage_data": {
            "summary": current_linkage_data["header"],
            "bin_data": current_linkage_data["bin_data"]
        }
    })


def modify_repeating_linkage(iterator, current_linkage_data):
    iterator["Linkage_data"] = {
        "summary": current_linkage_data["header"],
        "bin_data": current_linkage_data["bin_data"]
    }


def generate_json_output(binsDirectory):
    output = {
        "date_last_updated": date.today().strftime("%m/%d/%Y"),
        "database_name": "torsions_z_score_database",
        "data": None
    }
    output_data = []
    for currentFile in os.scandir(binsDirectory):
        sugar_info = parse_input_filename(currentFile.name)
        current_linkage_data = process_current_dotdata_file(currentFile.path)
        if current_linkage_data["header"][
                "count_mean"] is not None and current_linkage_data["header"][
                    "count_stdev"] is not None:
            match_donor_sugar = next((item for item in output_data
                                      if item["donor"] == sugar_info["donor"]),
                                     None)
            if match_donor_sugar is not None:
                list_of_acceptor_sugars = match_donor_sugar["acceptor"]
                match_acceptor_sugar = next(
                    (item for item in list_of_acceptor_sugars
                     if item["sugar"] == sugar_info["acceptor"]), None)
                if match_acceptor_sugar is not None:
                    list_of_linkages = match_acceptor_sugar["Linkage"]
                    matching_linkage = next(
                        (item for item in list_of_linkages
                         if item["donor_end"] == sugar_info["Linkage"]
                         ["donor_end"] and item["acceptor_end"] ==
                         sugar_info["Linkage"]["acceptor_end"]), None)
                    if matching_linkage is not None:
                        print(
                            "WARNING: WE SHOULD HAVE NOT GOTTEN HERE(LINE: 190)\nPrevious entry is going to be overwritten with the current one"
                        )
                        modify_repeating_linkage(matching_linkage,
                                                 current_linkage_data)
                    else:
                        append_output_for_unique_linkage(
                            match_acceptor_sugar, sugar_info["Linkage"],
                            current_linkage_data)
                else:
                    append_output_for_unique_acceptor_sugar(
                        match_donor_sugar, sugar_info["acceptor"],
                        sugar_info["Linkage"], current_linkage_data)
            else:
                output_data.append(
                    generate_output_for_unique_donor_sugar(
                        sugar_info, current_linkage_data))
        else:
            continue

    output["data"] = output_data

    return output


def main():
    scriptDirectory = get_path_of_directory_file_is_located_in()
    binsDirectory = os.path.join(scriptDirectory, input_bins_relative_path)
    exportJSON = generate_json_output(binsDirectory)

    outputDirectoryPath = os.path.join(scriptDirectory, output_directory)
    if os.path.exists(outputDirectoryPath):
        shutil.rmtree(outputDirectoryPath)
    CreateFolder(outputDirectoryPath)

    with open(
            os.path.join(outputDirectoryPath,
                         "privateer_torsions_z_score_database.json"),
            "w",
            encoding="utf-8",
    ) as export_json_file:
        json.dump(exportJSON, export_json_file, indent=3, ensure_ascii=False)


if __name__ == "__main__":
    main()