import os, csv
import matplotlib.pyplot as plt
import numpy as np


def get_path_of_directory_file_is_located_in():
    directory = os.path.dirname(os.path.abspath(__file__))
    return directory


def main():
    scriptDirectory = get_path_of_directory_file_is_located_in()

    bin = 2
    bin_size = int(360 / bin)

    torsion_data = {}
    total_no = 0
    filename = "NAG-NAG.csv"
    path = os.path.join(scriptDirectory, filename)

    with open(path, "r") as f:
        data = csv.reader(f)
        for index, entry_data in enumerate(data):
            if index == 0: continue
            linkage = entry_data[7]
            total_no += 1
            if linkage in torsion_data.keys():
                if entry_data[5] and entry_data[6]:
                    torsion_data[linkage]['phi'].append(float(entry_data[5]))
                    torsion_data[linkage]['psi'].append(float(entry_data[6]))
            else:
                torsion_data[linkage] = {'phi': [], 'psi': []}

    for linkage in torsion_data.keys():
        if linkage == "1-4":
            fig, ax = plt.subplots()
            hist, xbins, ybins, im = ax.hist2d(
                torsion_data[linkage]['phi'],
                torsion_data[linkage]['psi'],
                bins=(bin_size, bin_size),
                range=np.array([(-180, 180), (-180, 180)]),
                cmap=plt.get_cmap('gist_heat_r'))
            output = []
            count = []
            for i in range(len(ybins) - 1):
                for j in range(len(xbins) - 1):
                    current_item = [
                        xbins[j], xbins[j + 1], ybins[i], ybins[i + 1],
                        hist.T[i, j]
                    ]
                    # print(current_item)
                    output.append(current_item)
                    count.append(hist.T[i, j])

            sum_sq = 0
            sum_ = 0

            for x in count:
                sum_sq += x**2
                sum_ += x

            num = 0
            mean = sum_sq / sum_
            for x in count:
                num += x * ((x - mean)**2)

            stddev = num / (sum_ - 1)

            name = filename.split('.')[0]
            name_split = name.split('-')
            output_name = f"{name_split[0]}-{linkage}-{name_split[1]}"

            with open(os.path.join(scriptDirectory, f"{output_name}.data"),
                      'w') as output_file:
                output_file.write(f"{mean} {stddev}\n")
                for x in output:
                    for y in x:
                        output_file.write(str(y))
                        output_file.write(' ')
                    output_file.write('\n')


if __name__ == "__main__":
    main()