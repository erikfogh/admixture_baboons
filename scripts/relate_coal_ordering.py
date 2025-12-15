# Loading packages
import pandas as pd
import argparse
import tskit
import os

# Arg definition
parser = argparse.ArgumentParser(description="Analyze coal ordering from tree sequence")
parser.add_argument("-t", help="input tree")
parser.add_argument("-p", help="poplabel path")
# parser.add_argument("-i", help="focal ID")
parser.add_argument("-o", help="name and path to write the finished file")


def coalescence_ordering(ts, highest_ID, sample_counts, i_mapping):
    df_l = []
    c = 0
    for tree in ts.trees():
        if c % 2500 == 0:
            print(c)
        c += 1
        l_l = []
        l_c = []
        # Setup of base information
        for n in tree.timeasc()[:highest_ID]:
            l_l.append([n, tree.num_sites, tree.span, tree.interval[0]])
            l_c.append([])
        for n in tree.timeasc()[highest_ID:]:
            tree_time = tree.time(n)
            c_samples = (pd.Series([x for x in tree.samples(n)]))
            c_sample_counts = c_samples.map(i_mapping).value_counts()
            for p in c_sample_counts.index:
                if c_sample_counts[p] > sample_counts[p]/2:
                    for i in c_samples:
                        if p not in l_l[i]:
                            l_l[i].append(p)
                            l_c[i].append(tree_time)
        for n in tree.timeasc()[:highest_ID]:
            l_l[n].extend(l_c[n])
        df_l.append(pd.DataFrame(l_l))
    return pd.concat(df_l)


# # OLD implementation. Function for coal ordering
# def coalescence_ordering(tree, IDs, sample_counts):
#     df_list = []
#     for i in IDs:
#         pop_list = []
#         gen_list = []
#         current_node = i
#         while tree.depth(current_node) > 0:
#             # Find parent node
#             parent_node = tree.parent(current_node)
#             c_samples = pd.Series([x for x in tree.samples(current_node)])
#             c_sample_counts = c_samples.map(i_mapping).value_counts()
#             # If pop is not already added to list, add pop and note coal time when it surpasses 50 %
#             for p in c_sample_counts.index:
#                 if p  not in pop_list and c_sample_counts[p] > sample_counts[p]/2:
#                     pop_list.append(p)
#                     gen_list.append(tree.time(current_node))
#             current_node = parent_node
#         d = {"ID": i, "sites": tree.num_sites, "span": tree.span, "start": tree.interval[0]}
#         for i in range(len(pop_list)):
#             d["coal_{}".format(i)] = pop_list[i]
#         for i in range(len(gen_list)):
#             d["coal_date_{}".format(i)] = gen_list[i]
#         df_list.append(pd.DataFrame(d, index=[i]))
#     return pd.concat(df_list)


# Parsing args and constant defs
args = parser.parse_args()
treepath = args.t
poplabel_path = args.p
basepath = os.path.dirname(args.t)
treename = os.path.basename(args.t)


# Loading poplabels
poplabels = pd.read_csv(poplabel_path, sep=" ")

# Setup based on poplabels
if "hap" in treepath:
    sample_counts = poplabels["GROUP"].value_counts()
    i_mapping = {}
    for i, row in poplabels.iterrows():
        i_mapping[i] = row.GROUP
    highest_ID = poplabels.index.max()+1
else:
    sample_counts = poplabels["GROUP"].value_counts()*2
    i_mapping = {}
    for i, row in poplabels.iterrows():
        i_mapping[i*2] = row.GROUP
        i_mapping[i*2+1] = row.GROUP
    highest_ID = poplabels.index.max()*2+2


# Loading and running through the trees
ts = tskit.load(treepath)

print("Starting coalescence ordering with")
print("It has the following sample counts {}".format(sample_counts))
print("There are {} trees".format(len(ts.trees())))
cols = ["ID", "sites", "span", "start"]
for i in range(len(sample_counts)):
    cols.append("coal_{}".format(i))
for i in range(len(sample_counts)):
    cols.append("coal_date_{}".format(i))

df = coalescence_ordering(ts, highest_ID, sample_counts, i_mapping)
df.columns = cols

df.to_csv("{}/{}_coal_orders.txt".format(basepath, treename[:-6]), index=False)
