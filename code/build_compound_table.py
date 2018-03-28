import requests

def get_ftrd():
    """Get featuretablerawdata from database."""

    url = "http://coleykursarlab.chpc.utah.edu/api/featuretablerawdata/"
    ftrd = []
    count = 1

    while True:
        if count % 10 == 0:
            print(count, "pages processed")

        resp = requests.get(url)
        data = resp.json()
        ftrd += data["results"]

        if data["next"] is None:
            print(count, "pages total")
            break
        else:
            url = data["next"]

        count += 1

    return ftrd

def build_feature_table(raw_data, mz_error=0.01, rt_error=0.3):
    """Group features across samples.

    Arguments:
    raw_data -- feature table raw data downloaded from database
    mz_error -- maximum mass difference in Da (default 0.01)
    rt_error -- maximum retention time difference in minutes (default 0.3)
    """

    feature_table = {}
    current_feature_number = 0
    for datum in raw_data:
        pcid = datum["PC_ID"]
        mz = datum["MZ"]
        rt = datum["RT"]
        sample = datum["sample"]
        matched = False
        for feature in feature_table.keys():
            if abs(feature_table[feature]["avg_mz"] - mz) <= mz_error and \
                abs(feature_table[feature]["avg_rt"] - rt) <= rt_error:
                feature_table[feature]["mz"] += [mz]
                feature_table[feature]["rt"] += [rt]
                feature_table[feature]["PC_IDs"] += [pcid]
                feature_table[feature]["sample"] += [sample]
                feature_table[feature]["avg_mz"] = sum(feature_table[feature]["mz"])/ \
                 len(feature_table[feature]["mz"])
                feature_table[feature]["avg_rt"] = sum(feature_table[feature]["rt"])/ \
                 len(feature_table[feature]["rt"])
                feature_table[feature]["TIC"] += [datum["TIC"]]
                matched = True
                break
        if not matched:
            current_feature_number += 1
            feature_table[current_feature_number] = {
                "mz": [mz],
                "rt": [rt],
                "PC_IDs": [pcid],
                "sample": [sample],
                "TIC": [datum["TIC"]],
                "avg_mz": mz,
                "avg_rt": rt
                }

    return feature_table

def build_pc_id(feature_table):
    """Return dictionary of pcgroups and associated features.

    Arguments:
    feature_table -- output from build_feature_table function
    """

    pc_id_table = {}

    for feature in feature_table.keys():
        for pcid in set(feature_table[feature]["PC_IDs"]):
            tics = [feature_table[feature]["TIC"][i] for i, x in \
             enumerate(feature_table[feature]["PC_IDs"]) if x == pcid]
            if pcid in pc_id_table:
                pc_id_table[pcid]["features"] += [feature]
                pc_id_table[pcid]["TICs"] += [sum(tics) / len(tics)]
            else:
                pc_id_table[pcid] = {
                    "features": [feature],
                    "TICs": [sum(tics) / len(tics)]
                }

    return pc_id_table

def cosine_score(compound_1, compound_2, abundance, features):
    """Return cosine score of two vectors.

    Arguments:
    compound_1 -- dictionary entry containing named lists with feature and
     abundance values of components of first vector to be compared
    compound_2 -- dictionary entry containing named lists with feature and
     abundance values of components of second vector to be compared
    abundance -- name of list containing abundance values
    features -- name of list containing feature names
    """

    tics_1 = []
    tics_2 = []
    for feature in compound_1['features']:
        if feature in compound_2['features']:
            tics_1 += [float(compound_1[abundance][compound_1[features].index(feature)])]
            tics_2 += [float(compound_2[abundance][compound_2[features].index(feature)])]
    if len(tics_1) > 0:
        cosine_score = sum([tics_1[x] * tics_2[x] for x in range(len(tics_1))]) / \
         ((sum([x ** 2 for x in compound_1[abundance]]) ** (0.5)) * \
         (sum([x ** 2 for x in compound_2[abundance]]) ** (0.5)))
    else:
        cosine_score = 0

    return cosine_score

def build_compound_table(pc_id_table, min_cos_score=0.5):
    """Match pcgroups across samples based on shared features.

    Arguments:
    pc_id_table -- dictionary of pcgroups and associated features
     built with build_pc_id function
    min_cos_score -- pcgroups with cosine score above this cutoff
     are grouped as the same compound (default 0.5)
     """

    compound_table = {}
    sorted_pc_ids = sorted(pc_id_table.keys(), key=lambda k: len(pc_id_table[k]["features"]),\
     reverse=True)
    current_compound_number = 0

    for pcid in sorted_pc_ids:
        matched = False
        for compound in compound_table:
            if cosine_score(pc_id_table[pcid], compound_table[compound],\
             abundance="TICs", features="features") > min_cos_score:
                for feature in pc_id_table[pcid]["features"]:
                    if feature in compound_table[compound]["features"]:
                        compound_table[compound]["TICs"]\
                         [compound_table[compound]["features"].index(feature)] += \
                         pc_id_table[pcid]["TICs"][pc_id_table[pcid]["features"].index(feature)]
                    else:
                        compound_table[compound]["features"] += [feature]
                        compound_table[compound]["TICs"] +=\
                         [pc_id_table[pcid]["TICs"][pc_id_table[pcid]["features"].index(feature)]]
                matched = True
                break
        if not matched:
            current_compound_number += 1
            compound_table[current_compound_number] = {
                "features": list(pc_id_table[pcid]["features"]),
                "TICs": list(pc_id_table[pcid]["TICs"])
            }
    for compound in compound_table:
        total_tic = sum(compound_table[compound]["TICs"])
        compound_table[compound]["feature_pcts"] = [x / total_tic \
        for x in compound_table[compound]["TICs"]]
        max_tic = max(compound_table[compound]["feature_pcts"])
        compound_table[compound]["rel_abund"] = [x / max_tic \
        for x in compound_table[compound]["feature_pcts"]]

    return compound_table

def fill_compounds(filled_features, compound_table):
    """Return dictionary containing compounds found in each sample.

    Arguments:
    filled_features -- file created using R code as in fill_peaks_all_inga.R containing
     abundance values for all features found in each sample.
    compound_table -- dictionary of compounds and associated features created by
     build_compound_table function.
    """

    filled_compounds = {}

    for sample in filled_features:
        for compound in compound_table:
            top_feature = compound_table[compound]["features"][compound_table[compound]["rel_abund"].index(1.0)]
            if top_feature in filled_features[sample]["feature_number"]:
                actual_rt = filled_features[sample]["actual_rt"]\
                [filled_features[sample]["feature_number"].index(top_feature)]
            else:
                continue
            if not 1.0 < actual_rt < 32.0:
                continue
            shared_features = [x for i, x in enumerate(filled_features[sample]["feature_number"]) if \
             x in compound_table[compound]["features"] and \
             abs(filled_features[sample]["actual_rt"][i] - actual_rt) < 0.1]
            major_features = [x for i, x in enumerate(compound_table[compound]["features"]) if \
             compound_table[compound]["rel_abund"][i] > 0.8]
            shared_major_features = [x for x in shared_features if \
             x in major_features]
            if len(shared_major_features) == len(major_features):
                sample_compound = {
                    "features": shared_features,
                    "TICs": [filled_features[sample]["TIC"]\
                     [filled_features[sample]["feature_number"].index(x)] for x in shared_features],
                    "actual_rt": [filled_features[sample]["actual_rt"]\
                     [filled_features[sample]["feature_number"].index(x)] for x in shared_features]
                }
                top_feature_tic = sample_compound["TICs"][sample_compound["features"].index(top_feature)]
                sample_compound["rel_abund"] = [x / top_feature_tic for x in sample_compound["TICs"]]
                actual_rel_abund = [compound_table[compound]["rel_abund"][compound_table[compound]["features"].index(x)]\
                 for x in sample_compound["features"]]
                feature_within_range = [0.33 < actual_rel_abund[i] / sample_compound["rel_abund"][i] < 3 for i, x in enumerate(actual_rel_abund)]
                if False in feature_within_range:
                    sample_compound_2 = {
                        "features": [x for i, x in enumerate(sample_compound["features"]) if feature_within_range[i]],
                         "TICs": [x for i, x in enumerate(sample_compound["TICs"]) if feature_within_range[i]]
                    }
                    sample_compound = sample_compound_2

                #sample_compound["rel_abund"] = [x / max(sample_compound["TICs"]) for x in sample_compound["TICs"]]
                #shared_rel_abund = \
                # [sample_compound["rel_abund"][sample_compound["features"].index(x)] for x in shared_major_features]
                
                #if False in [x > 0.33 for x in shared_rel_abund]:
                #    low_abund_feature = shared_major_features[shared_rel_abund.index(min(shared_rel_abund))]
                #    low_abund = sample_compound["TICs"][sample_compound["features"].index(low_abund_feature)]
                #    sample_compound_2 = {
                #        "features": [x for i, x in enumerate(sample_compound["features"]) if \
                #         sample_compound["TICs"][i] < 3 * low_abund],
                #         "TICs": [x for i, x in enumerate(sample_compound["TICs"]) if \
                #         sample_compound["TICs"][i] < 3 * low_abund]
                #    }
                    
                if cosine_score(compound_table[compound], sample_compound, abundance="TICs", features="features")\
                >= 0.3:
                    if sample in filled_compounds:
                        filled_compounds[sample]["compound"] += [compound]
                        filled_compounds[sample]["TIC"] += [sum(sample_compound["TICs"])]
                    else:
                        filled_compounds[sample] = {
                            "compound": [compound],
                            "TIC": [sum(sample_compound["TICs"])]
                        }

                    fields = ["TIC", "actual_rt", "feature_number"]

                    used_feature_idx = [filled_features[sample]["feature_number"].index(k) for k in sample_compound["features"]]
                    for field in fields:
                        filled_features[sample][field] = [x for i, x in \
                        enumerate(filled_features[sample][field]) if i not in used_feature_idx]

    return filled_compounds
