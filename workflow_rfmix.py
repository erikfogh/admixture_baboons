from gwf import Workflow, AnonymousTarget
import os
from groups import Group
import pandas as pd

gwf = Workflow()

### Input paths
# Uses the baboondiversity environment. Requirements are mostly bcftools and rfmix.

path_to_vcfs = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/chr{}/chr{}.phased.rehead.vcf.gz"
sim_base = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/sim_haptools/"
path_to_sims = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/sim_haptools/"
x_all_path = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/chrX_with_males/chrX_diploid_all_nomiss.vcf.gz"
genetic_map = "/home/eriks/baboondiversity/data/PG_panu3_recombination_map/mikumi_pyrho_genetic_map_chr{}.txt"
path_to_output = "steps/rfmix_gen100/"
base_path = os.getcwd()
autosomes = list(range(1, 21)) #+ ["X"]

meta_data_samples = pd.read_csv("data/Papio_metadata_with_clustering_sci.txt", sep=" ")

# Creating lists of the various inputs needed to run prep and rfmix.
# For the first iteration, query is everyone except the references and gelada
# but this could be altered.

full_list = ['Cynocephalus, Central Tanzania', 'Anubis, Kenya', 'Kindae, Zambia',
             'Hamadryas, Ethiopia', 'Anubis, Tanzania',
             'Cynocephalus, Western Tanzania', 'Papio, Senegal', 'Ursinus, Zambia',
             'Anubis, Ethiopia']

ref_name_list = [["tanzania_focus", ['Ursinus, Zambia', 'Kindae, Zambia',
            'Hamadryas, Ethiopia', 'Papio, Senegal']],
                ["eth_olive_focus", ['Hamadryas, Ethiopia', 'Papio, Senegal',
            'Cynocephalus, Central Tanzania', 'Anubis, Tanzania']],
                ["tanzania_gm_focus", ['Gog Woreda, Gambella region, Ethiopia', 'Mikumi, Tanzania']],
                ["tanzania_pure_focus", ['Gog Woreda, Gambella region, Ethiopia',
                                         'Mahale, Tanzania', 'Mikumi, Tanzania']]]


map_inputs = []

for n in ref_name_list:
    os.makedirs(path_to_output+"/"+n[0], exist_ok=True)
    os.makedirs(path_to_output+"/"+n[0]+"_plain", exist_ok=True)
    if n[0] == "tanzania_gm_focus" or n[0] == "tanzania_pure_focus":
        meta_data_samples_sub = meta_data_samples.loc[meta_data_samples.Origin.isin(n[1])]
        query_samples = meta_data_samples.loc[~(meta_data_samples.Origin.isin(n[1])) &
                                            (meta_data_samples.Origin != "Gelada, Captive")]
        pop_df = pd.DataFrame({"PDGP_ID": meta_data_samples_sub.PGDP_ID,
                               "population": meta_data_samples_sub.Origin})
        pop_df.to_csv(path_to_output+"/"+n[0]+"/ref_names.txt",
                      index=False, header=False, sep="\t")
    
        ref_samples_f = meta_data_samples.loc[~(meta_data_samples.Origin.isin(n[1])) &
                                            (meta_data_samples.Origin != "Gelada, Captive") &
                                            (meta_data_samples.Sex == "F")]
        pop_df = pd.DataFrame({"PDGP_ID": meta_data_samples_sub.PGDP_ID,
                               "population": meta_data_samples_sub.Origin})
        pop_df.to_csv(path_to_output+"/"+n[0]+"/female_ref_names.txt",
                      index=False, header=False, sep="\t")
    else:
        meta_data_samples_sub = meta_data_samples.loc[meta_data_samples.C_origin.isin(n[1])]
        query_samples = meta_data_samples.loc[~(meta_data_samples.C_origin.isin(n[1])) &
                                            (meta_data_samples.C_origin != "Gelada, Captive")]
        pop_df = pd.DataFrame({"PDGP_ID": meta_data_samples_sub.PGDP_ID,
                               "population": meta_data_samples_sub.C_origin})
        pop_df.to_csv(path_to_output+"/"+n[0]+"/ref_names.txt",
                      index=False, header=False, sep="\t")
    
        ref_samples_f = meta_data_samples.loc[~(meta_data_samples.C_origin.isin(n[1])) &
                                            (meta_data_samples.C_origin != "Gelada, Captive") &
                                            (meta_data_samples.Sex == "F")]
        pop_df = pd.DataFrame({"PDGP_ID": meta_data_samples_sub.PGDP_ID,
                               "population": meta_data_samples_sub.C_origin})
        pop_df.to_csv(path_to_output+"/"+n[0]+"/female_ref_names.txt",
                      index=False, header=False, sep="\t")
    d = {}
    d["run_name"] = n[0]
    d["ref_samples"] =list(meta_data_samples_sub.PGDP_ID)
    d["query_samples"] = list(query_samples.PGDP_ID)
    map_inputs.append(d)
    
    # # Hap approach - obsolete, left if I need it for other workflows.
    # x_hap_IDs, x_hap_pops = [], []
    # for i, row in meta_data_samples.iterrows():
    #     if row.Sex == "F":
    #         x_hap_IDs.append(row.PGDP_ID + "_1")
    #         x_hap_IDs.append(row.PGDP_ID + "_2")
    #         x_hap_pops.append(row.C_origin), x_hap_pops.append(row.C_origin)
    #     else:
    #         x_hap_IDs.append(row.PGDP_ID)
    #         x_hap_pops.append(row.C_origin)

    # hap_df = pd.DataFrame({"PDGP_ID": x_hap_IDs,
    #                        "population": x_hap_pops})
    # ref_df = hap_df.loc[hap_df.population.isin(n[1])]
    # query_df = hap_df.loc[~hap_df.population.isin(n[1])]
    # hap_df.to_csv(path_to_output+"/"+n[0]+"/hap_ref_names.txt",
    #               index=False, header=False, sep="\t")
    # d = {}
    # d["run_name"] = n[0]
    # d["ref_samples"] =list(ref_df.PDGP_ID)
    # d["query_samples"] = list(query_df.PDGP_ID)
    # hap_inputs.append(d)


### Gwf functions

def prep_rfmix(run_name, ref_samples, query_samples, out_suffix, chr_list, path_to_output):
    output_ref = path_to_output + run_name + "/" + out_suffix+"_ref.bcf"
    output_query = path_to_output + run_name + "/" + out_suffix+"_query.bcf"
    input_match = path_to_vcfs.format("*", chr_list)
    if out_suffix == "X_female":
        chr_iter = "X"
        inputs = path_to_vcfs.format(chr_iter, chr_iter)
    elif out_suffix == "X_all":
        inputs = x_all_path
        input_match = x_all_path
    else:
        s_chr = chr_list.split("..")
        chr_iter = list(range(int(s_chr[0][1:]), int(s_chr[1][:-1])))
        inputs = [path_to_vcfs.format(x, x) for x in chr_iter]
    ref = ",".join(ref_samples)
    query = ",".join(query_samples)
    outputs = [output_ref, output_query]
    options = {
        "cores": 2,
        "memory": "30g",
        "walltime": "4:00:00",
        "account": "baboondiversity"
    }
    spec = """
    bcftools concat {input_match} | bcftools view  -s {ref}  -q 0.01:minor \
    -O b -o {output_ref} --force-samples
    bcftools index {output_ref}
    bcftools concat {input_match} | bcftools view -s {query} -q 0.01:minor \
    -O b -o {output_query} --force-samples
    bcftools index {output_query}
    """.format(input_match=input_match, ref=ref, query=query,
               output_ref=output_ref, output_query=output_query)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def prep_rfmix_sim(run_name, path_to_sims, out_suffix, chr_list, path_to_output):
    output_query = path_to_output + out_suffix+"_query.bcf"
    inputs = [path_to_sims]
    outputs = [output_query]
    options = {
        "cores": 2,
        "memory": "30g",
        "walltime": "4:00:00",
        "account": "baboondiversity"
    }
    spec = """
    bcftools concat {path_to_sims} | bcftools view -q 0.01:minor \
    -O b -o {output_query} --force-samples
    bcftools index {output_query}
    """.format(path_to_sims=path_to_sims,
            output_query=output_query)
    # print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def rfmix(chrom, query, reference, sample_map, genetic_map, output_path, gen=100):
    # Run with "low" memory, then with higher memory to complete all jobs
    output = output_path + "chr" + str(chrom)
    m = genetic_map
    inputs = [query, reference, m]
    outputs = [output + ".msp.tsv"]
    options = {
        "cores": 10,
        "memory": "200g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    rfmix -f {} -r {} -m {} -g {} -o {} --chromosome=chr{} -e 3 -G {} --reanalyze-reference
    """.format(query, reference, sample_map, m, output, chrom, gen)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def rfmix_plain(chrom, query, reference, sample_map, genetic_map, output_path, gen=100):
    # Run with "low" memory, then with higher memory to complete all jobs
    output = output_path + "chr" + str(chrom)
    m = genetic_map
    inputs = [query, reference, m]
    outputs = [output + ".msp.tsv"]
    options = {
        "cores": 10,
        "memory": "200g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    rfmix -f {} -r {} -m {} -g {} -o {} --chromosome=chr{} -e 5 -G {}
    """.format(query, reference, sample_map, m, output, chrom, gen)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

### Function calls


if not os.path.exists(path_to_output + "aut_genetic_map.txt"):
    print(path_to_output + "aut_genetic_map.txt")
    df_l = []
    for a in range(1, 21):
        recomb_df = pd.read_csv(genetic_map.format(a), sep=" ")
        df_l.append(recomb_df)
    (pd.concat(df_l)[["chromosome", "position", "Genetic_Map(cM)"]]).to_csv(path_to_output + "aut_genetic_map.txt",
                             sep=" ", index=False)
    print("Created recomb df for aut")

gwf.map(prep_rfmix, map_inputs, name="autosomes",
        extra={"out_suffix": "aut", "chr_list": "{1..20}",
               "path_to_output": path_to_output})

# Autosomal loop

for i in range(len(ref_name_list)):
    n = ref_name_list[i][0]
    file_name = path_to_output + n
    gwf.map(rfmix, autosomes, name="rfmix_"+n,
                extra={"query": file_name+"/aut_query.bcf", "reference": file_name+"/aut_ref.bcf",
                       "sample_map": file_name+"/ref_names.txt",
                       "genetic_map": path_to_output + "aut_genetic_map.txt",
                       "output_path": file_name+"/"})
    gwf.map(rfmix_plain, autosomes, name="rfmix_plain_"+n,
                extra={"query": file_name+"/aut_query.bcf", "reference": file_name+"/aut_ref.bcf",
                       "sample_map": file_name+"/ref_names.txt",
                       "genetic_map": path_to_output + "aut_genetic_map.txt",
                       "output_path": file_name+"_plain/"})


if not os.path.exists(path_to_output + "X_genetic_map.txt"):
    print(path_to_output + "X_genetic_map.txt")
    df_l = []
    for a in ["X"]:
        recomb_df = pd.read_csv(genetic_map.format(a), sep=" ")
        df_l.append(recomb_df)
    (pd.concat(df_l)[["chromosome", "position", "Genetic_Map(cM)"]]).to_csv(path_to_output + "X_genetic_map.txt",
                             sep=" ", index=False)
    print("Created recomb df for chrX")

gwf.map(prep_rfmix, map_inputs, name="chrX_female",
        extra={"out_suffix": "X_female", "chr_list": ["X"],
               "path_to_output": path_to_output})

gwf.map(prep_rfmix, map_inputs, name="chrX_all",
        extra={"out_suffix": "X_all", "chr_list": "X",
               "path_to_output": path_to_output})


for i in range(len(ref_name_list)):
    n = ref_name_list[i][0]
    file_name = path_to_output + n
    gwf.map(rfmix, ["X"], name="rfmix_X_female"+n,
                extra={"query": file_name+"/X_female_query.bcf", "reference": file_name+"/X_female_ref.bcf",
                       "sample_map": file_name+"/female_ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": file_name+"/female_"})
    gwf.map(rfmix, ["X"], name="rfmix_X_all"+n,
                extra={"query": file_name+"/X_all_query.bcf", "reference": file_name+"/X_all_ref.bcf",
                       "sample_map": file_name+"/ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": file_name+"/all_"})
    gwf.map(rfmix, ["X"], name="rfmix_X_female_ref"+n,
                extra={"query": file_name+"/X_all_query.bcf", "reference": file_name+"/X_female_ref.bcf",
                       "sample_map": file_name+"/female_ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": file_name+"/female_ref_"})
    # Plain runs
    gwf.map(rfmix, ["X"], name="rfmix_X_female_plain"+n,
                extra={"query": file_name+"/X_female_query.bcf", "reference": file_name+"/X_female_ref.bcf",
                       "sample_map": file_name+"/female_ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": file_name+"_plain/female_"})
    gwf.map(rfmix, ["X"], name="rfmix_X_all_plain"+n,
                extra={"query": file_name+"/X_all_query.bcf", "reference": file_name+"/X_all_ref.bcf",
                       "sample_map": file_name+"/ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": file_name+"_plain/all_"})
    gwf.map(rfmix, ["X"], name="rfmix_X_female_ref_plain"+n,
                extra={"query": file_name+"/X_all_query.bcf", "reference": file_name+"/X_female_ref.bcf",
                       "sample_map": file_name+"/female_ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": file_name+"_plain/female_ref_"})
    

# # Simulation runs. It is meant to replicate the Tanzania analysis, so it uses the same reference as that.

os.makedirs(path_to_output+"sim_results/", exist_ok=True)
os.makedirs(path_to_output+"sim_results_plain/", exist_ok=True)


sim_settings = [{"run_name": "aut_50", "model": "gog_mikumi_50gen.dat", "chrom": "8"},
                  {"run_name": "X_50", "model": "gog_mikumi_50gen.dat", "chrom": "X"},
                  {"run_name": "X_250", "model": "gog_mikumi_250gen.dat", "chrom": "X"}]

sim_list_dict = []
for s in sim_settings:
    d = {}
    d["run_name"] = s["run_name"]
    d["path_to_sims"] = path_to_sims+s["run_name"]+s["model"].split(".")[0]+".vcf.gz"
    d["out_suffix"] = "chr{}_sim_{}".format(s["chrom"], s["model"].split(".")[0])
    d["chr_list"] = "{"+s["chrom"]+".."+s["chrom"]+"}"
    sim_list_dict.append(d)

gwf.map(prep_rfmix_sim, sim_list_dict, name="sims",
        extra={"path_to_output": path_to_output+"sim_results/"})

gen_list = [50, 100, 250, 500, 1000]

for g in gen_list:
    output_add = "infer{}_".format(g)
    gwf.map(rfmix, ["8"], name="rfmix_aut_sim_gog_mik"+output_add,
                extra={"query": path_to_output+"sim_results/chr8_sim_gog_mikumi_50gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/aut_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/ref_names.txt",
                       "genetic_map": path_to_output + "aut_genetic_map.txt",
                       "output_path": path_to_output+"sim_results/"+output_add,
                       "gen": g})

    gwf.map(rfmix, ["X"], name="rfmix_X50_sim_gog_mik"+output_add,
                extra={"query": path_to_output+"sim_results/chrX_sim_gog_mikumi_50gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/X_all_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": path_to_output+"sim_results/gen50_"+output_add,
                       "gen": g})

    gwf.map(rfmix, ["X"], name="rfmix_X50_sim_gog_mik_female"+output_add,
                extra={"query": path_to_output+"sim_results/chrX_sim_gog_mikumi_50gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/X_female_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/female_ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": path_to_output+"sim_results/gen50_female_"+output_add,
                       "gen": g})


    gwf.map(rfmix, ["X"], name="rfmix_X250_sim_gog_mik"+output_add,
                extra={"query": path_to_output+"sim_results/chrX_sim_gog_mikumi_250gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/X_all_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": path_to_output+"sim_results/gen250_"+output_add,
                       "gen": g})

    gwf.map(rfmix, ["X"], name="rfmix_X250_sim_gog_mik_female"+output_add,
                extra={"query": path_to_output+"sim_results/chrX_sim_gog_mikumi_250gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/X_female_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/female_ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": path_to_output+"sim_results/gen250_female_"+output_add,
                       "gen": g})

    # Plain runs

    gwf.map(rfmix_plain, ["8"], name="rfmix_aut_sim_gog_mik_plain"+output_add,
                extra={"query": path_to_output+"sim_results/chr8_sim_gog_mikumi_50gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/aut_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/ref_names.txt",
                       "genetic_map": path_to_output + "aut_genetic_map.txt",
                       "output_path": path_to_output+"sim_results_plain/"+output_add,
                       "gen": g})

    gwf.map(rfmix_plain, ["X"], name="rfmix_X50_sim_gog_mik_plain"+output_add,
                extra={"query": path_to_output+"sim_results/chrX_sim_gog_mikumi_50gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/X_all_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": path_to_output+"sim_results_plain/gen50_"+output_add,
                       "gen": g})

    gwf.map(rfmix_plain, ["X"], name="rfmix_X50_sim_gog_mik_female_plain"+output_add,
                extra={"query": path_to_output+"sim_results/chrX_sim_gog_mikumi_50gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/X_female_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/female_ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": path_to_output+"sim_results_plain/gen50_female_"+output_add,
                       "gen": g})


    gwf.map(rfmix_plain, ["X"], name="rfmix_X250_sim_gog_mik_plain"+output_add,
                extra={"query": path_to_output+"sim_results/chrX_sim_gog_mikumi_250gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/X_all_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": path_to_output+"sim_results_plain/gen250_"+output_add,
                       "gen": g})

    gwf.map(rfmix_plain, ["X"], name="rfmix_X250_sim_gog_mik_female_plain"+output_add,
                extra={"query": path_to_output+"sim_results/chrX_sim_gog_mikumi_250gen_query.bcf",
                       "reference": path_to_output+"tanzania_focus/X_female_ref.bcf",
                       "sample_map": path_to_output+"tanzania_focus/female_ref_names.txt",
                       "genetic_map": path_to_output + "X_genetic_map.txt",
                       "output_path": path_to_output+"sim_results_plain/gen250_female_"+output_add,
                       "gen": g})
