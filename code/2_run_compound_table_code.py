import imp
ct = imp.load_source('ct', './/code//build_compound_table.py')

####### load peak table created in xcms_peak_picking.R #######
with open(".//results//all_peaks_long.csv", 'r') as load_peaks:
    peaks = []
    for line in load_peaks:
        line_1 = line.replace('"','').replace('\n','').split(',')
        if line_1[0] == "sample":
            continue
        peaks.append({
            "RT": float(line_1[2]),
            "PC_ID": line_1[4],
            "TIC": float(line_1[5]),
            "species_code_sample": line_1[1],
            "MZ": float(line_1[3]),
            "sample": line_1[0]
        })

####### match features across samples using build_feature_table function #######
# each unique feature is assigned a unique number.
# set mz tolerance (in Da) and rt tolerance (in minutes)
feature_table = ct.build_feature_table(peaks, mz_error=0.01, rt_error=0.3)
# write csv of feature table
with open(".//results//features_matched_across_samples.csv", 'w') as file1:
    file1.write("feature_number,sample,rt,mz,PC_ID,TIC\n")
    for feature in feature_table:
        for sample in feature_table[feature]["sample"]:
            file1.write("%d,%s,%f,%f,%s,%f\n" % (feature, sample,\
             feature_table[feature]["rt"][feature_table[feature]["sample"].index(sample)],\
             feature_table[feature]["mz"][feature_table[feature]["sample"].index(sample)],\
             feature_table[feature]["PC_IDs"][feature_table[feature]["sample"].index(sample)],\
             feature_table[feature]["TIC"][feature_table[feature]["sample"].index(sample)]))

####### Build table detailing which features are associated with each PC_ID
pc_id_table = ct.build_pc_id(feature_table)

####### Match PC_IDs across samples and assign unique compound numbers.
compound_table = ct.build_compound_table(pc_id_table, min_cos_score=0.5)
# write csv of compound table
with open(".//results//compound_feature_table.csv", "w") as file1:
    file1.write("compound_number,feature_number,TIC,mz,rt\n")
    for compound in compound_table:
        for feature in compound_table[compound]["features"]:
            file1.write("%d,%d,%f,%f,%f\n" % (compound, feature, \
            compound_table[compound]["TICs"][compound_table[compound]["features"].index(feature)], \
            feature_table[feature]["avg_mz"], feature_table[feature]["avg_rt"]))