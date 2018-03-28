import imp
ct = imp.load_source('ct', './/code//build_compound_table.py')

####### load filled features created in 'fill_peaks.R':
with open(".//results//filled_features.csv", "r") as file1:
    filled_features = {}
    for row in file1:
        sample = row.split(",")[4].split("\n")[0].replace('"', '')
        if sample == 'sample_name':
            continue
        else:
            if sample in filled_features:
                if int(row.split(",")[0]) not in filled_features[sample]["feature_number"]:
                    filled_features[sample]["feature_number"] += [int(row.split(",")[0])]
                    filled_features[sample]["TIC"] += [float(row.split(",")[1])]
                    filled_features[sample]["actual_rt"] += [float(row.split(",")[3])]
            else:
                filled_features[sample] = {
                    "feature_number": [int(row.split(",")[0])],
                    "TIC": [float(row.split(",")[1])],
                    "actual_rt": [float(row.split(",")[3])]
                }

####### Use filled features and fill_compounds function to identify which compounds are present in each sample
filled_comps = ct.fill_compounds(filled_features, compound_table)

# write csv of final compound table
with open(".//results//filled_compounds.csv", "w") as file1:
    file1.write("species_sample,compound_number,TIC\n")
    for sample in filled_comps:
        for i, compound in enumerate(filled_comps[sample]["compound"]):
            file1.write("%s,%d,%f\n" % (sample, compound, \
            filled_comps[sample]["TIC"][i]))