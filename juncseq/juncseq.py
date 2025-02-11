# juncseq.py
import argparse
import pyranges as pr
import pandas as pd
import numpy as np
import glob
import os


def read_list(list_file):
    """Read a list from a file and return it as a set."""
    with open(list_file, 'r') as f:
        result = {line.strip() for line in f if line.strip()}
    return result


def subset_gtf(gtf_file, gene_list):
    """
    Subset a GTF file based on a given gene list.

    Args:
        gtf_file (str): Path to the input GTF file.
        gene_list (set): Set of gene names to subset by.
    """
    annotations = pr.read_gtf(gtf_file)

    filtered_annotations = annotations[annotations.transcript_id.isin(gene_list['transcript_id'])]
    filtered_annotations = filtered_annotations[filtered_annotations.Feature == "exon"]
    filtered_columns = filtered_annotations[ ['gene_name', 'transcript_id', 'exon_number']]
    return filtered_columns

def get_gene_coordinates(df, buffer):
    """
    Compute the lowest start value and highest end value for each gene in the dataset.

    Args:
        df (pd.DataFrame): A DataFrame containing columns 'Chromosome', 'Start', 'End', 'Strand', and 'gene_name'.

    Returns:
        pd.DataFrame: A DataFrame with one row per gene, including the chromosome, strand, 
                      lowest start, and highest end values.
    """
    # Group by 'gene_name' and compute required values
    result = (
        df.groupby("gene_name")
        .agg(
            Chromosome=("Chromosome", "first"),  # Assuming one chromosome per gene
            Strand=("Strand", "first"),          # Assuming one strand per gene
            Start=("Start", "min"),
            End=("End", "max"),
            Exons=("exon_number", "max"),
        )
        .reset_index()
    )

    result["Start"] = result["Start"] - buffer
    result["End"] = result["End"] + buffer
    return result

def extract_junctions(junction_files, exon_coordinates):
    result_df = pd.DataFrame()
    files_processed = 0
    gene_coordinates = get_gene_coordinates(exon_coordinates, 200)

    for file in junction_files:
        # Read the junction file and add a header
        all_junctions = pd.read_csv(file, sep='\t', header=None)
        headers = ['chromosome', 'intronStart', 'intronEnd', 'strand', 'canonicalMotif', 'annotated', 'uniqueReads', 'multimappedReads', 'overhang'] 
        all_junctions.columns = headers

        for index, gene in gene_coordinates.iterrows():
            in_region = filter_junctions(all_junctions, gene)
            in_region = annotate_junctions(in_region, exon_coordinates[exon_coordinates['gene_name'] == gene['gene_name']])
            in_region['JAF_start'], in_region['JAF_end'], in_region['totalReads_start'], in_region['totalReads_end'] = calculate_JAF(in_region)

            # Filter out low frequency junctions
            high_JAF = in_region[(in_region['JAF_start'] >= 0.05) | (in_region['JAF_end'] >= 0.05)]
            
            # Make columns int for reporting
            columns_to_convert = ["start_exonNumber", "end_exonNumber", "start_distance", "end_distance"]
            high_JAF[columns_to_convert] = high_JAF[columns_to_convert].astype(int)
            
            # Annotate junctions with splicing explanation
            high_JAF['classification'] = classify_splicing(high_JAF)

            # Add a column with the filename (without the path and extension)
            high_JAF["Filename"] = os.path.basename(file).replace(".hs38.SJ.out.tab", "")

            # Append to the result DataFrame
            result_df = pd.concat([result_df, high_JAF])

        # Reporting log
        files_processed += 1
        if files_processed % 10 == 0:
            print(f"Processed {files_processed}/{len(junction_files)} junction files")


    return result_df

def filter_junctions(junctions, gene):
    if gene['Strand'] == '+':
        gene_strand = 1
    else:
        gene_strand = 2

    # Apply filtering based on the chromosome and position ranges
    filtered = junctions[
        (junctions['uniqueReads'] != 0) &
        (junctions['overhang'] >= 10) &
        (junctions['strand'] == gene_strand) &
        (junctions['chromosome'] == gene['Chromosome']) &
        (
            ((junctions['intronStart'] >= gene['Start']) & (junctions['intronStart'] <= gene['End'])) |  # Intron start position in region 
            ((junctions['intronEnd'] >= gene['Start']) & (junctions['intronEnd'] <= gene['End']))    # OR intron end position in region
        )
    ]
    return filtered

def calculate_JAF(junctions):
    # Don't allow for the start and end junctions to be the same (for JAF)
    for index, row in junctions.iterrows():
        if row['start_exonNumber'] == row['end_exonNumber']:
            junctions.at[index, 'end_exonNumber'] += 1
    
    exonStart_agg = (
        junctions.groupby("start_exonNumber")
        .agg(
            totalReads_start=("uniqueReads", "sum")
        )
        .reset_index()
    )

    exonEnd_agg = (
        junctions.groupby("end_exonNumber")
        .agg(
            totalReads_end=("uniqueReads", "sum")
        )
        .reset_index()
    )

    junctions = junctions.merge(exonStart_agg, on='start_exonNumber', how='left')
    junctions = junctions.merge(exonEnd_agg, on='end_exonNumber', how='left')     
    JAF_start = round(junctions["uniqueReads"]/junctions["totalReads_start"],2)
    JAF_end = round(junctions["uniqueReads"]/junctions["totalReads_end"],2)

    return JAF_start, JAF_end, junctions['totalReads_start'], junctions['totalReads_end']

def annotate_junctions(junctions, exons): 
    exons['exon_number'] = exons['exon_number'].astype(int)   
    if exons['Strand'].values[0] == "-":
        # Reverse exon_number labelling 
        n_exons = exons['exon_number'].max()
        exons['exon_number'] = n_exons + 1 - exons['exon_number']

    intron_start = exons[['End', 'exon_number']]
    intron_start['End'] += 1 # Shift to line up with junction value
    intron_end = exons[['Start', 'exon_number']]

    # Rename and merge files to annotate junctions with exon number
    intron_start.columns = ['intronStart', 'start_exonNumber']
    intron_end.columns = ['intronEnd', 'end_exonNumber']

    merged = junctions.merge(intron_start, on='intronStart', how='left')
    merged = merged.merge(intron_end, on='intronEnd', how='left')   
    
    # Apply the function only to rows where `start_exonNumber` is NaN
    na_mask_start = merged['start_exonNumber'].isna()
    start_distance, start_closestExonNumber = find_closest_exon(
        merged.loc[na_mask_start, 'intronStart'].to_numpy(),        # Junction starts
        intron_start['intronStart'].to_numpy(),        # Intron starts
        intron_start['start_exonNumber'].to_numpy(),   # Exon numbers
        "start"
    )

    # Apply the function only to rows where `end_exonNumber` is NaN
    na_mask_end = merged['end_exonNumber'].isna()
    end_distance, end_closestExonNumber = find_closest_exon(
        merged.loc[na_mask_end, 'intronEnd'].to_numpy(),        # Junction starts
        intron_end['intronEnd'].to_numpy(),        # Intron starts
        intron_end['end_exonNumber'].to_numpy(),   # Exon numbers
        "end"
    )

    # Update the `merged` DataFrame
    merged.loc[na_mask_start, 'start_exonNumber'] = start_closestExonNumber
    merged.loc[na_mask_end, 'end_exonNumber'] = end_closestExonNumber

    if exons['Strand'].values[0] == "+":
        merged.loc[na_mask_start, 'start_distance'] = start_distance
        merged.loc[na_mask_end, 'end_distance'] = end_distance
    else:
        merged.loc[na_mask_start, 'start_distance'] = -start_distance
        merged.loc[na_mask_end, 'end_distance'] = -end_distance

    merged.fillna(0, inplace=True) # Write distances as 0 if it is a canonical junction

    return merged


def find_closest_exon(junction_coord, intron_coord, exon_numbers, region):
    """
    Calculate the closest `exonNumber` and distance to junction
    """
    if region == "start":
        diff = intron_coord - junction_coord[:, None] 
    else:
        diff = junction_coord[:, None] - intron_coord 

    closest_indices = np.abs(diff).argmin(axis=1)
    signed_diffs = diff[np.arange(len(junction_coord)), closest_indices]
    closest_exons = exon_numbers[closest_indices]
    return signed_diffs, closest_exons


def classify_splicing(junctions):
    classification = []
    for index, junction in junctions.iterrows():
        # Both junctions annotated
        if junction['start_distance'] == 0 and junction['end_distance'] == 0: 
            # Sequential exons
            exon_span = junction['end_exonNumber'] - junction['start_exonNumber']
            if exon_span == 1:
                classification.append('Normal')
            else:
                classification.append('Skipping (' + str(exon_span-1) + ')')
        
        # Start junction annotated
        elif junction['start_distance'] == 0:
            if junction['end_distance'] > 0:
                classification.append('3prime_exonExtension')
            else: 
                classification.append('3prime_exonTruncation')

        # End junction annotated
        elif junction['end_distance'] == 0:
            if junction['start_distance'] > 0:
                classification.append('5prime_exonExtension')
            else: 
                classification.append('5prime_exonTruncation')

        else:        
            classification.append('Else')

    return classification

def main():
    parser = argparse.ArgumentParser(description="JuncSeq: Analyze splice junctions for aberrant splicing.")
    parser.add_argument("--patient-list", required=True, help="Path to the file containing the list of junction files.")
    parser.add_argument("--gtf", required=True, help="Path to the gene annotation.")
    parser.add_argument("--gene-list", required=False, help="Path to the file containing the list of genes.")
    parser.add_argument("--vcf", required=False, help="Path to the VCF file containing variants.")
    parser.add_argument("--output", required=True, help="Path to the output directory.")

    args = parser.parse_args()

    # Your processing logic here
    print(f"Analyzing gene list: {args.gene_list}")
    genes = pd.read_csv(args.gene_list)

    print(f"Using gene annotation defined by: {args.gtf}")
    gene_gtf = subset_gtf(args.gtf, genes)
    
    junction_files = read_list(args.patient_list)
    print(f"Analysing junctions for {len(junction_files)} files defined in: {args.patient_list}")
    all_junctions = extract_junctions(junction_files, gene_gtf.df) # Make GTF a df
    all_junctions.to_csv(f"{args.output}/all_junctions.tsv", sep='\t', index=False)

    filtered_junctions = all_junctions[all_junctions["classification"].str.contains("Skipping", na=False)]
    filtered_junctions.to_csv(f"{args.output}/skipping.tsv", sep='\t', index=False)

if __name__ == "__main__":
    main()