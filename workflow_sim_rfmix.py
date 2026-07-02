from gwf import Workflow, AnonymousTarget
import os
from groups import Group

gwf = Workflow()

### Input paths
# Uses the baboondiversity environment. Requirements are mostly bcftools and rfmix.

path_to_vcfs = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/chr{}/chr{}.phased.rehead.vcf.gz"
sim_base = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/sim_haptools/"
path_to_gog_mik = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/sim_haptools/chr{}_gog_mik.vcf.gz"
#path_to_mik_gog = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/sim_haptools/chr{}_mik_gog.vcf.gz"
x_all_path = "/home/eriks/baboondiversity/data/PG_panu3_phased_chromosomes_4_7_2021/chrX_with_males/chrX_diploid_all_nomiss.vcf.gz"
genetic_map = "/home/eriks/baboondiversity/data/PG_panu3_recombination_map/mikumi_pyrho_genetic_map_chr{}.txt"
path_to_output = "steps/rfmix_gen100/"
base_path = os.getcwd()
autosomes = list(range(1, 21)) #+ ["X"]


### Gwf functions

def haptools_sim(model, map_file, chr_number, in_vcf, sample_tsv, out_vcf):
    inputs=[model, map_file, in_vcf, sample_tsv]
    outputs = [out_vcf]
    options = {
        "cores": 10,
        "memory": "120g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = """
    haptools simgenotype --model {model}  --mapdir {map_file} --chroms {chr_number} \
    --ref_vcf {in_vcf} --sample_info {sample_tsv} --out {out_vcf}
    """.format(model=model, map_file=map_file, chr_number=chr_number,
               in_vcf=in_vcf, sample_tsv=sample_tsv, out_vcf=out_vcf)
    print(spec)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


### Function calls

# # Simulation runs. It is meant to replicate the Tanzania analysis, so it uses the same reference as that.

input_settings = [{"run_name": "aut_50", "model": "gog_mikumi_50gen.dat", "chrom": 8},
                  {"run_name": "X_50", "model": "gog_mikumi_50gen.dat", "chrom": "X"},
                  {"run_name": "X_250", "model": "gog_mikumi_250gen.dat", "chrom": "X"}]

sim_list_dict = []
for s in input_settings:
    d = {}
    d["model"] = sim_base+s["model"]
    d["map_file"] = sim_base+"map_files/"
    d["chr_number"] = s["chrom"]
    if s["chrom"] == "X":
        d["in_vcf"] = x_all_path
    else:
        d["in_vcf"] = path_to_vcfs.format(s["chrom"], s["chrom"])
    d["sample_tsv"] = sim_base+"aut.tsv"
    d["out_vcf"] = sim_base+s["run_name"]+s["model"].split(".")[0]+".vcf.gz"
    sim_list_dict.append(d)

print(sim_list_dict)

gwf.map(haptools_sim, sim_list_dict)

# gwf.map(prep_rfmix_sim, [{"run_name": "aut_sim_mik_gog", "path_to_sims": path_to_mik_gog}], name="autosome_sim_mik_gog",
#         extra={"out_suffix": "aut_sim_50gen_mik_gog", "chr_list": "{8..8}",
#                "path_to_output": path_to_output})
# gwf.map(rfmix, ["8"], name="aut_sim_mik_gog",
#                 extra={"query": path_to_output+"aut_sim_mik_gog/aut_sim_50gen_mik_gog_query.bcf",
#                        "reference": path_to_output+"tanzania_focus/aut_ref.bcf",
#                        "sample_map": path_to_output+"tanzania_focus/ref_names.txt",
#                        "genetic_map": path_to_output + "aut_genetic_map.txt",
#                        "output_path": path_to_output+"aut_sim_mik_gog/"})
