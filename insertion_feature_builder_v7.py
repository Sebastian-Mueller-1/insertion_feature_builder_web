import re

import numpy as np
import pandas as pd
from pybedtools import BedTool


def find_overlap_genomic_feature_identifier(df_to_build, source_df):
    """This function will add the genomic identifier of interest to the first column of the output df. This might be Ds_GFP alleles in the context of the project this code was written for. This id column could be some other genomic feature of interest. But it assumes the feature of interest is the name listed in the the ID='' information in attributes field"""

    # extract the feature id, in this case the Ds GFP allele name
    df_to_build["feature_of_interest_ID"] = source_df["source_attributes"].map(
        lambda x: x.split("=")[1].split(";")[0]
    )

    # get start and end
    df_to_build["feature_of_interest_ID_start"] = source_df["source_start"]
    df_to_build["feature_of_interest_ID_end"] = source_df["source_end"]

    return df_to_build


def only_canonical(overlap_df):
    """This function will filter the overlap df to only contain overlaps that are associated with the cannonical transctipt id, this is not necessary for the other function to run, but will produce conflated column entries that mix and match non canonical and canonical transcript overlaps. For example calling the genomic region classifer might output yes for 3' UTR even if overlap is not on 3 prime UTR of canonical transcript but on the 3 prime UTR of some other transcript"""

    # make temp df to hold canonical transctip IDs
    temp_df = pd.DataFrame(columns=["Canonical_transcript_ID"])

    # use lambda function to pull out the transcript ids for canonical only transcripts
    temp_df["Canonical_transcript_ID"] = overlap_df["ref_attributes"].apply(
        lambda x: x.split("ID=")[1].split(";")[0]
        if "canonical_transcript=1" in x
        else None
    )

    # need to drop non rows that the apply function makes when nothing is found
    temp_df = temp_df.dropna()
    temp_df_list = temp_df[
        "Canonical_transcript_ID"
    ].tolist()  # produce a list to use in next step regex

    temp_df_set = set(temp_df_list)

    # use pandas regex str contains and join method to make a long 'or' regex to check each row of the overlap df attributes column to check whether that overlap is canonical associated or not
    canonical_only_overlap_df = overlap_df[
        overlap_df["ref_attributes"]
        .apply(lambda x: re.search("Z[^;]*_T[^;]*", x))
        .apply(lambda match: match.group(0) if match else None)
        .isin(temp_df_set)
    ]

    # canonical_only_overlap_df = overlap_df[overlap_df['ref_attributes'].str.contains("|".join(temp_df_list))]

    # make a copy to avoid pandas weird df slice error
    canonical_only_overlap_df = canonical_only_overlap_df.copy()

    # add the canonical transcript IDs as a column in the overlap df, this will be needed for using the find CDS length and Distance from CDS start function
    canonical_only_overlap_df["Canonical_transcript_ID"] = canonical_only_overlap_df[
        "ref_attributes"
    ].apply(
        lambda x: re.search("Z[^;]*_T[^;]*", x).group(0)
        if re.search("Z[^;]*_T[^;]*", x)
        else None
    )

    # add the gene id from these canonical transcript id by stripping the transcript tag
    canonical_only_overlap_df["Gene_ID"] = canonical_only_overlap_df[
        "Canonical_transcript_ID"
    ].map(lambda x: x.split("_")[0])

    return canonical_only_overlap_df


def find_overlap_gene_info(df_to_build, overlap_df, find_length):
    '''This function finds genes in the overlap df and adds them to the df_to_build using a pd merge method. It will map many to many for genomic features of interest that may overlap with multiple genes. It will enter these as additional columns it assumses youve already produced the output dataframe with genomic features of interest. This function also accepts an boolean argument "find_length, if True it will calculate the length of the gene overlaps."'''

    # subset the overlap df for relevant info for this function and apply a filter to only keep overlaps with gene feature ID. This is done to make the pd meerge call simpler
    overlap_df_subset = overlap_df[
        [
            "feature_of_interest_ID",
            "ref_feature",
            "ref_attributes",
            "ref_start",
            "ref_end",
            "ref_seqname",
            "ref_strand",
            "source_strand",
        ]
    ]
    overlap_df_subset = overlap_df_subset[overlap_df_subset["ref_feature"] == "gene"]

    # perform merge and allow for many to many key matches
    df_to_build = pd.merge(
        df_to_build, overlap_df_subset, on="feature_of_interest_ID", how="inner"
    )

    # perform a maping function to pull a string out
    df_to_build["Gene_ID"] = df_to_build["ref_attributes"].map(
        lambda x: x.split("=")[1].split(";")[0]
    )

    # if find length True map out length
    if find_length is True:
        df_to_build["Gene_length"] = df_to_build.apply(
            lambda row: row["ref_end"] - row["ref_start"] + 1, axis=1
        )

    # remove junk columns used for the merge and string extraction and rename one columns
    df_to_build = df_to_build.drop(
        ["ref_attributes", "ref_feature", "ref_end", "ref_start"], axis=1
    )
    df_to_build = df_to_build.rename(
        columns={
            "ref_seqname": "feature_of_interest_ID_chromosome",
            "ref_strand": "Gene_ID_strand",
            "source_strand": "feature_of_interest_ID_strand",
        }
    )

    return df_to_build


def check_feature_conditions_dooner_ds_gfp_specific(row, feature):
    """check the specific 8bp orientation condition for dooner ds gfp. + + or - + you need to make sure end of ds gfp insertion 8 bp is in span of feature. + - or - - you need to make sure start of ds gfp allele"""

    # Skip rows where 'feature' is already False
    if row[feature] is False:
        return row

    # check conditions from doc string
    if (row["feature_of_interest_ID_strand"] == "+" and row["ref_strand"] == "+") or (
        row["feature_of_interest_ID_strand"] == "-" and row["ref_strand"] == "+"
    ):
        if not (
            row["ref_start"] <= row["feature_of_interest_ID_end"] <= row["ref_end"]
        ):
            row[feature] = False
            # print('Dooner orientation expeception for:', row)

    elif (row["feature_of_interest_ID_strand"] == "+" and row["ref_strand"] == "-") or (
        row["feature_of_interest_ID_strand"] == "-" and row["ref_strand"] == "-"
    ):
        if not (
            row["ref_start"] <= row["feature_of_interest_ID_start"] <= row["ref_end"]
        ):
            row[feature] = False
            # print('Dooner orientation expeception for:', row)

    return row


def feature_class_overlap(df_to_build, canonical_only_overlap_df, feature_array):
    """This function will produce true false columns if genomic of feature of interest overlaps with genomic feature in ref annotation gff3 file. feature_array column takes a list of the string names of genomic feature (e.g. "exon","CDS", etc). Assumes overlap df is canonical only"""

    canonical_only_overlap_df = (
        canonical_only_overlap_df.copy()
    )  # use this copy feature to avoid pandas df slicing copy issues

    # make a Gene_ID column for the the overlap df needed for this feature overlap function
    canonical_only_overlap_df["Gene_ID"] = canonical_only_overlap_df[
        "ref_attributes"
    ].map(lambda x: x.split("=")[1].split("_")[0])

    for feature in feature_array:
        # subset the overlap df to contain only relevant columns
        overlap_df_subset = canonical_only_overlap_df[
            [
                "Gene_ID",
                "ref_feature",
                "feature_of_interest_ID",
                "ref_strand",
                "ref_start",
                "ref_end",
            ]
        ]

        # filter the subset to contain only the feature of interest (e.g. exon, cds, UTR, etc.)
        overlap_df_subset = overlap_df_subset[
            overlap_df_subset["ref_feature"] == feature
        ]

        # use merge to get genes that have insertions on feature of interest
        df_to_build = pd.merge(
            df_to_build,
            overlap_df_subset,
            on=["Gene_ID", "feature_of_interest_ID"],
            how="outer",
        )

        # built the nicely named feature column with True or False if insertion is present in overlap
        df_to_build[feature] = (
            df_to_build["ref_feature"]
            .str.contains(feature)
            .astype(str)
            .replace("nan", False)
        )

        # Apply the check_feature_conditions function
        df_to_build = df_to_build.apply(
            lambda row: check_feature_conditions_dooner_ds_gfp_specific(row, feature),
            axis=1,
        )

        # drop the unnecessary columns that were made in merge
        df_to_build.drop(
            ["ref_feature", "ref_start", "ref_end", "ref_strand"], axis=1, inplace=True
        )

    return df_to_build


def find_canonical_transcript_id(df_to_build, canonical_only_overlap_df):
    """Find the canonical transcript id of each genetic feature of interest ID. assumes a canonical only_overlap df and gene ID"""

    # make a filtered subset df of just gene id and transcript id
    filtered_mRNA_only_df = canonical_only_overlap_df[
        canonical_only_overlap_df["ref_feature"] == "mRNA"
    ]
    subset_df = filtered_mRNA_only_df[
        ["Gene_ID", "Canonical_transcript_ID", "feature_of_interest_ID"]
    ]

    # run merge on two dfs using gene id and feature of interest ID as key matches. This should hand cases where there may be repeat genes or for the repeat feature of interest id
    df_to_build = pd.merge(
        df_to_build, subset_df, on=["Gene_ID", "feature_of_interest_ID"], how="left"
    )

    return df_to_build


def find_cds_metrics(df_to_build, canonical_annotation_df):
    """This function will find the CDS length for each canonincal transcrip Gene in the def_to_build, assumes canonical only overlap df. requires a canonical transcrip id in df_to_build. run the find canonical transcript id function requires feature class overlap to be run."""

    def cds_length_and_distance_from_start_calculator(row, cds_only_bed):
        """This is a helper function that takes a canonical transcrip id and finds the actual length for that cds and the distance from start for genomic feature of interest"""

        # use try clause since there will be some NaN cells and we want those filled with None
        try:
            canonical_transcript_id = row["Canonical_transcript_ID"]
            feautre_of_interest_start = row["feature_of_interest_ID_start"]
            feature_of_interest_end = row["feature_of_interest_ID_end"]

            # filter canonical cds bed to only be rows with a given transctip id
            filtered_cds_only_bed = cds_only_bed.filter(
                lambda x: canonical_transcript_id in x[8]
            ).saveas()

            # grab strand the cds is running on. will need this for getting orientation start right for Ds_gfp start
            strand = filtered_cds_only_bed[0].strand

            # run merge wich collapses the bed to just contain start and finish and chr columns
            merge_cds = filtered_cds_only_bed.merge()

            # init running sum variable to handle cases where no valid start cds overlap and will be used to find distance from start calculation in valid casses
            running_sum = 0
            bp_distance_from_start = 0

            # need to cases for neg and pos strand because start orientation is opposite. Ds gfps have repeating 8bp seq preceeding and flanking actual insertion. so there is a logic behind if it inserts at start or end of Ds deppending on orientation. Recommend drawing a diagram to understand why positive cds use end and negative cds use start

            if strand == "+":  # positive CDS strand
                # iteration across each row in merge cds bed object
                for feature in merge_cds:
                    # run conditional to to see if ds gfp end is in the range (start to finish) of the current merge cds feature
                    if int(feature_of_interest_end) in range(
                        int(feature[1]), int(feature[2])
                    ):
                        # if it is calculate bp distance from start by getting distance from current feature start and adding it to the running sum
                        bp_distance_from_start = (
                            int(feature_of_interest_end) - int(feature[1]) + running_sum
                        )

                    # calculate running sum using the bed too feature length method
                    running_sum += feature.length

            else:  # negative CDS strand case
                # iterate backwards across the merge CDS object will follow above logic but use the start of the ds gfp instead of start

                for feature in reversed(merge_cds):
                    if int(feautre_of_interest_start) in range(
                        int(feature[1]), int(feature[2])
                    ):
                        bp_distance_from_start = (
                            int(feature[2])
                            - int(feautre_of_interest_start)
                            + running_sum
                        )
                    running_sum += feature.length

            # add this conditional to handle cases where genomic feature does not overlap with CDS. This will make the bp related columns in fine df state false to conform to naming standards in this program

            if row["CDS"] == "True":
                return [
                    running_sum,
                    bp_distance_from_start,
                    bp_distance_from_start / running_sum,
                ]
            return [running_sum, None, None]
        # except handles issues with the strand = filtered_cds_only_bed[0].strand assignment when the filter returns none for a no canonical overlap genomic feature of interest
        except IndexError:
            return [None, None, None]
        except TypeError:
            return [None, None, None]

    # make canonical only bed object and filter for CDS only rows. Do this to use convienient Bed methods for bioinformatics. Drop two extra columns so that the canonical annotation df will follow needed form for bed tool.

    canonical_only_df_for_bed = canonical_annotation_df.copy()
    canonical_only_df_for_bed = canonical_only_df_for_bed.drop(
        ["ref_attributes", "Gene_ID"], axis=1
    )

    # make bed object and filter for CDS only
    # canonical_only_bed = BedTool.from_dataframe(canonical_only_df_for_bed)
    canonical_only_bed = BedTool.from_dataframe(canonical_only_df_for_bed)
    cds_only_bed = canonical_only_bed.filter(lambda x: x[2] == "CDS").saveas()

    # run the helper function for each transcript
    df_to_build["processed_CDS_info"] = df_to_build.apply(
        lambda row: cds_length_and_distance_from_start_calculator(row, cds_only_bed),
        axis=1,
    )

    # expand the processed cds column to three different columns
    df_to_build[
        [
            "CDS_length",
            "bp_distance_from_start_of_CDS",
            "bp_proportion_from_start_of_CDS",
        ]
    ] = df_to_build["processed_CDS_info"].apply(pd.Series)

    df_to_build = df_to_build.drop("processed_CDS_info", axis=1)

    df_to_build.replace(np.nan, "False", inplace=True)

    return df_to_build


def fast_feature_length_finder(canonical_only_annotation_df, feature):
    single_feature_df = canonical_only_annotation_df.copy()
    single_feature_df = single_feature_df.drop(["ref_attributes", "Gene_ID"], axis=1)
    single_feature_df = single_feature_df[single_feature_df["ref_feature"] == feature]
    print(single_feature_df)

    single_feature_df.sort_values(
        by=["Canonical_transcript_ID", "ref_start"], inplace=True
    )

    # Calculate the CDS length for each row
    single_feature_df["CDS_Length"] = (
        single_feature_df["ref_end"] - single_feature_df["ref_start"] + 1
    )

    # Group by Transcript_ID and aggregate the needed information
    single_feature_df = (
        single_feature_df.groupby("Canonical_transcript_ID")
        .agg(
            start_position=("ref_start", "min"),
            end_position=("ref_end", "max"),
            cds_length=("CDS_Length", "sum"),
            strand=(
                "ref_strand",
                "first",
            ),  # Assuming strand is the same for all rows of a transcript
        )
        .reset_index()
    )

    # Set Transcript_ID as the index if needed
    single_feature_df.set_index("Canonical_transcript_ID", inplace=True)

    return single_feature_df


def multiple_associated_genes_warning(df_to_build):
    """This function will look at current df and output a warning file that has all the genomic features of interest that are associated with multiple genes"""

    def get_genes(group):
        """This helper actually finds the repeats and values of gene associate ID"""

        if len(group) >= 2:
            genes = group["Gene_ID"].tolist()
            return pd.Series([genes])
        else:
            return pd.Series([None])

    grouped_genes = df_to_build.groupby("feature_of_interest_ID", as_index=False).apply(
        get_genes
    )
    grouped_genes = grouped_genes.dropna().reset_index(drop=True)

    grouped_genes.columns = ["Feature_of_interest_ID", "Gene_IDs"]

    return grouped_genes


def remove_specified_repeat_rows(remove_df, df_to_build):
    """This function removes rows with specified feature of interest and gene id that is passed initially to the program. If the file is empty it will return the original df unchanged"""

    for index, row in remove_df.iterrows():
        # strip spaces
        feature_id = str(row.iloc[0]).strip()
        gene_ID = str(row.iloc[1]).strip()

        # build a match condition boolean mask for all the rows
        match_condition = (df_to_build["feature_of_interest_ID"] == feature_id) & (
            df_to_build["Gene_ID"] == gene_ID
        )

        # use ~ to negative the bool mask to remove the rows that have postive match
        df_to_build = df_to_build[~match_condition]

    return df_to_build


def lift_off(lift_off_df, df_to_build):
    """This function will do a lift off for the current working df"""

    # rename lift off df to match naming of df to build
    old_name1 = lift_off_df.columns[0]
    old_name2 = lift_off_df.columns[1]

    lift_off_df = lift_off_df.rename(
        columns={old_name2: "Gene_ID", old_name1: "Gene_ID_v3"}
    )

    # want to collapse the df down incase a gene_id has multiple gene IDs in other version
    collapsed_lift_off_df = (
        lift_off_df.groupby("Gene_ID")["Gene_ID_v3"].agg(list).reset_index()
    )

    # perform merge to get new version of gene ID
    df_to_build = df_to_build.merge(
        collapsed_lift_off_df[["Gene_ID", "Gene_ID_v3"]], on="Gene_ID", how="left"
    )

    # make NaN False for posterity
    df_to_build["Gene_ID_v3"] = df_to_build["Gene_ID_v3"].fillna("False")

    return df_to_build


def process_data(
    annotation_path, insertion_path, remove_csv, lift_off_txt
):
    # need this df to build some genomic feature of interest information for down stream functions. In context of this project it is the Ds gfp allele gff file
    source_df = pd.read_csv(
        insertion_path,
        sep="\t",
        comment="#",
        names=[
            "source_seqname",
            "source_source",
            "source_feature",
            "source_start",
            "source_end",
            "source_score",
            "source_strand",
            "source_phase",
            "source_attributes",
        ],
    )

    # need this annotation df for down stream for cds metrics function call
    annotation_df = pd.read_csv(
        annotation_path,
        sep="\t",
        comment="#",
        names=[
            "ref_seqname",
            "ref_source",
            "ref_feature",
            "ref_start",
            "ref_end",
            "reg_score",
            "ref_strand",
            "ref_phase",
            "ref_attributes",
        ],
    )

    # make the remove feature df
    remove_df = pd.read_csv(
        remove_csv, names=["feature_ID", "Gene_ID", "Comment"], comment="#"
    )

    # lift_off_df = pd.read_csv(lift_off_txt, comment="#", sep="\t")

    # print(source_df)
    # print("--------------------------------------")
    # print(annotation_df)
    # print("--------------------------------------")
    # print(remove_df)
    # print("--------------------------------------")

    # ================ run bedtools overlap and make dfs ========================================================
    bed1 = BedTool(annotation_path)
    bed2 = BedTool(insertion_path)

    # Find overlapping regions
    overlap_regions = bed1.intersect(bed2, c=False, wo=True)
    # print(overlap_regions)

    # Initalize column names for overlap regions pandas df. Ref will be the annotation, and source will be your source gff3 you provide to find overlaps with annotation gff3
    column_names = [
        "ref_seqname",
        "ref_source",
        "ref_feature",
        "ref_start",
        "ref_end",
        "reg_score",
        "ref_strand",
        "ref_phase",
        "ref_attributes",
        "source_seqname",
        "source_source",
        "source_feature",
        "source_start",
        "source_end",
        "source_score",
        "source_strand",
        "source_phase",
        "source_attributes",
        "bed_tools_count",
    ]

    overlap_df = overlap_regions.to_dataframe(names=column_names)

    # Drop a special bed tools count that is added when calling the overlap method
    overlap_df = overlap_df.drop(["bed_tools_count"], axis=1)

    # This will add the genomic feature of interest ID from the attribute column of the source df to the master overlap output df. For example, it will take the ds gfp allele name and add it as a column for every overlap row. This method will be required to run for other functions to work

    overlap_df["feature_of_interest_ID"] = overlap_df["source_attributes"].map(
        lambda x: x.split("=")[1].split(";")[0]
    )

    # initalize output df. This will be the master df that gets all featurs added to it
    output_df = pd.DataFrame()

    output_df = find_overlap_genomic_feature_identifier(output_df, source_df)
    output_df = find_overlap_gene_info(output_df, overlap_df, True)
    overlap_df_canonical_only = only_canonical(overlap_df)
    output_df = feature_class_overlap(
        output_df,
        overlap_df_canonical_only,
        ["exon", "CDS", "five_prime_UTR", "three_prime_UTR"],
    )
    output_df = find_canonical_transcript_id(output_df, overlap_df_canonical_only)
    canonical_only_annotation_df = only_canonical(annotation_df)
    output_df = find_cds_metrics(output_df, canonical_only_annotation_df)
    multiple_warnings_df = multiple_associated_genes_warning(output_df)
    output_df = remove_specified_repeat_rows(remove_df, output_df)
    # output_df = lift_off(lift_off_df, output_df)

    # Return the processed DataFrames
    return output_df, multiple_warnings_df, remove_df


# if __name__ == "__main__":
#     # Just for direct script testing, paths are hardcoded
#     output_df, multiple_warnings_df, remove_df = process_data(
#         "/path/to/annotation.gff3",
#         "/path/to/insertion.gff3",
#         "/path/to/template.csv",
#         "/path/to/remove.csv",
#         "/path/to/liftover.txt"
#     )
#     # Example code to save these DataFrames to CSV files or further process them
