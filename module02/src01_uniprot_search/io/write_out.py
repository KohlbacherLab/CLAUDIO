
def write_sqcs_output(data, intra_only, filename, output_directory):
    # write output sqcs-file
    #
    # input data: pd.DataFrame, intra_only: bool, filename: str, output_directory: str
    # no return

    columns = ["pep_a", "pep_b", "pos_a", "pos_b", "seq_a", "seq_b",
               "unip_id_a", "unip_id_b", "gene_a", "gene_b", "pub", "XL_type"]
    output_path = f"{output_directory}{filename}.sqcs"
    if not intra_only:
        data[columns].to_csv(output_path, index=True)
    else:
        data = data[data["XL_type"] == "intra"]
        data[columns].to_csv(output_path, index=True)

        # Write issue output file containing all entries for which peptides are not directly contained in the full
        # sequence, for these a revision of the given/found uniprot ids might be useful, for example search for isoforms
        out_str = "unip_id,site_id,line,site_pos,pep,full_seq\n"
        c = 0
        uc = [0, []]
        for i, row in data.iterrows():
            if row["pep_a"] not in row["seq_a"]:
                out_str += f"{row['unip_id_a']},a,{i + 2},{row['pos_a']},{row['pep_a']},{row['seq_a']}\n"
                c += 1
                if row['pep_a'] not in uc[1]:
                    uc[0] += 1
                    uc[1].append(row['pep_a'])
            if row["pep_b"] not in row["seq_b"]:
                out_str += f"{row['unip_id_b']},b,{i + 2},{row['pos_b']},{row['pep_b']},{row['seq_b']}\n"
                c += 1
                if row['pep_b'] not in uc[1]:
                    uc[0] += 1
                    uc[1].append(row['pep_b'])
        out_str += f"=====================\nTotal count: {c}\nTotal unique peptide count: {uc[0]}"
        with open(f"{output_directory}{filename}_issue.csv", 'w') as f:
            f.write(out_str)