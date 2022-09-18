import os
from privateer import privateer_core as pvt

PDB_dir = '/home/harold/Dev/privateer_python/tests/test_data/5fjj.pdb'
# PDB_dir = '/home/harold/Dev/quick_python_utility_scripts/4ln8.pdb'
current_PDB = '4ln8'
try:
    glycosylation = pvt.GlycosylationComposition_memsafe(PDB_dir)
except RuntimeError:
    print(f"{current_PDB} caused a runtime error!")
except:
    print(os.error)

number_of_glycans = glycosylation.get_number_of_glycan_chains_detected()
OfflineTorsionsZScoreDB = pvt.OfflineTorsionsZScoreDatabase()

for glycanNo in range(number_of_glycans):

    glycan = glycosylation.get_glycan(glycanNo)
    numsugars = glycan.get_total_number_of_sugars()
    glycan_tmp = {}

    glycan_torsion_zscore = glycan.get_torsions_zscore_summary(
        OfflineTorsionsZScoreDB)
    print(glycan_torsion_zscore)
    # if (len(glycan_torsion_zscore)):
    #     print(glycan_torsion_zscore[0]["zscore"] is None)

    # print(glycan.get_total_of_glycosidic_bonds() == len(glycan_torsion_zscore))
    for sugar_index in range(numsugars):
        sugar = glycan.get_monosaccharide(numsugars - sugar_index - 1)
        linkage = sugar.get_sugar_linkage_info()
        glycan_summary = glycan.get_glycan_summary()  #
        sugar_name = sugar.get_name_short() + "-" + glycan_summary['RootInfo'][
            'RootSugarChainID'] + "-" + sugar.get_sugar_pdb_id().strip()
        sugar_id = sugar.get_sugar_id()
        glycan_tmp[sugar_id] = sugar_name

        for link in linkage:
            pass
            # print(glycan_tmp[link['connectedToSugarID']], glycan_tmp[sugar_id],
            #       link['linkageTorsions'])
            # print(link)
            # linkage_list.append({
            #     "Sugar 1": glycan_tmp[sugar_id],
            #     "Sugar 2": glycan_tmp[link['connectedToSugarID']],
            #     "Torsions": link['linkageTorsions']
            # })

    glycan_tmp = {}