import requests
import pandas as pd
import json
import argparse

def read_genes(file_path):
    with open(file_path, 'r') as file:
        genes = file.read().strip().split("\n")
    return genes

def query_enrichr(gene_list, background_list, description,gene_set_library):
    """
    Possible gene_set_libraries:
        - WikiPathways_2023_Human
        - Chromosome_Location
        - DisGeNET
        - Disease_Perturbations_from_GEO_down, Disease_Perturbations_from_GEO_upp
        - GO_Biological_Process_2023
        - GO_Cellular_Component_2023
        - GO_Molecular_Function_2023
        - GTEx_Tissue_Expression_Up/Down
        - GWAS_Catalog_2023
        - KEGG_2021_Human
        - MSigDB_Hallmark_2020
        - OMIM_Expanded
    #We want wones focused on phenotypes, if possible
    """
    user_response = uploadGeneList(gene_list, description)
    background_response= uploadBackgroundGeneList(background_list)
    enrichment_response = testEnrichment(user_response,background_response, gene_set_library)
    #print(enrichment_response)

    """
    Code donated by chatGPT
        enrichr_url = "https://maayanlab.cloud/Enrichr/addList"
    gene_str = "\n".join(gene_list)
    background_str = "\n".join(background_list)
    
    # Posting gene list
    response = requests.post(enrichr_url, files={'list': (None, gene_str)})
    if not response.ok:
        raise Exception("Error in submitting gene list to Enrichr")
    result = response.json()
    user_list_id = result['userListId']
    
    # Enrichment analysis
    enrichr_url = f"https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType=genomic"
    response = requests.get(enrichr_url)
    if not response.ok:
        raise Exception("Error in performing enrichment analysis")
    result = response.json()
    """

    return enrichment_response


def uploadBackgroundGeneList(bg_list):
    base_url = "https://maayanlab.cloud/speedrichr"
    res = requests.post(
        base_url+'/api/addbackground',
        data=dict(background='\n'.join(bg_list)),
    )

    if res.ok:
        background_response = res.json()
        return background_response
    else:
        print("Error uploading background list")
        return None



def testEnrichment(userlist_response,background_response, gene_set_library):

    ## Get enrichment results with a background set:
    base_url = "https://maayanlab.cloud/speedrichr"

    res = requests.post(
            base_url+'/api/backgroundenrich',
            data=dict(
            userListId=userlist_response["userListId"],
            backgroundid=background_response["backgroundid"],
            backgroundType=gene_set_library,
            )
        )
    if res.ok:
        results = res.json()
        return results
    else:
         raise Exception("Error in performing enrichment analysis")

    #Get the enrichment results:
    """
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    gene_set_library = 'KEGG_2015'
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)
    print(data)
    """



def uploadGeneList(gene_list, description):
    """
    Uploads  gene list to test for enrichment
    """
    #Upload a gene list to analyze
    base_url = "https://maayanlab.cloud/speedrichr"
    res = requests.post(
        base_url+'/api/addList',
        files=dict(
        list=(None, '\n'.join(gene_list)),
        description=(None, description)
        ))

    if res.ok:
        userlist_response = res.json()
        return userlist_response
    else:
        raise Exception('Error uploading gene list')


def save_to_dataframe(result, output_path):
    df = pd.json_normalize(result)
    df.to_csv(output_path, index=False)


import pandas as pd

def parse_enrichr_results(enrichr_json,output_path,gene_set_lib):
    # Initialize an empty list to store the parsed data
    parsed_data = []

    # Iterate over each category in the JSON results
    for category, results in enrichr_json.items():
        for result in results:
            rank = result[0]
            term_name = result[1]
            p_value = result[2]
            odds_ratio = result[3]
            combined_score = result[4]
            adjusted_p_value = result[6]
            genes = ",".join(result[5])
            # Append the extracted data as a row
            parsed_data.append({
                "Rank": rank,
                "Term name": term_name,
                "P-value": p_value,
                "Odds ratio": odds_ratio,
                "Combined score": combined_score,
                "Adjusted p-value": adjusted_p_value,
                "Genes":genes
            })

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(parsed_data)
    df["gene_library"]=gene_set_lib
    #print(df)
    df.to_csv(output_path, index=False)
    return df

def main(genes_to_test_path, background_genes_path, output_path, gene_set_lib, uniq):
    genes_to_test = read_genes(genes_to_test_path)
    background_genes = read_genes(background_genes_path)
    if uniq:
      genes_to_test= list(set(genes_to_test))
      background_genes= list(set(background_genes))
    result = query_enrichr(genes_to_test, background_genes, "test",gene_set_lib)
    parse_enrichr_results(result, output_path,gene_set_lib)
    #save_to_dataframe(result, output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Query Enrichr for gene enrichment analysis")
    parser.add_argument('-t', '--test', required=True, help="Path to the file containing the list of genes to test")
    parser.add_argument('-b', '--background', required=True, help="Path to the file containing the background gene list")
    parser.add_argument('-o', '--output', required=True, help="Path to save the output CSV file")
    parser.add_argument('-g', '--gene_set_lib', default="ChEA_2022", help="Specify which gene set library to test- see the enrichr webiste for the names")
    parser.add_argument('-u', '--unique_only', default=False,action="store_true", help="Specify this if you are only wnating to test on unique sets of genes")
    args = parser.parse_args()

    main(args.test, args.background, args.output, args.gene_set_lib, args.unique_only)

