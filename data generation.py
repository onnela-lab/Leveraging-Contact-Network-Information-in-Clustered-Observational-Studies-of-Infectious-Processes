import os, sys, networkx as nx, numpy as np, pandas as pd, scipy as sy, copy
import scipy.io
import os.path
#os.chdir("C:/Users/Patrick/Desktop")
data_file = "data/datav4.0/datav4.0/"
network_files = "Data/1. Network Data/"
ID_files      = "Data/2. Demographics and Outcomes/"

df  = pd.read_stata(data_file+ID_files+"individual_characteristics.dta")

villages = list(range(1,78))
villages.remove(13)
villages.remove(22)

outcomes, As, PIDs, HHIDs, leaders = dict(), dict(), dict(), dict(), dict()
for village in villages:
    As[village] = np.loadtxt(data_file+network_files+"Adjacency Matrices/adj_allVillageRelationships_vilno_"+str(village)+".csv", delimiter=",")
    PIDs[village] = np.loadtxt(data_file+network_files+"Adjacency Matrix Keys/key_vilno_"+str(village)+".csv", dtype=int)
    HHIDs[village] = np.loadtxt(data_file+network_files+"Adjacency Matrix Keys/key_HH_vilno_"+str(village)+".csv", dtype=int)
    mf_filename = data_file+"Matlab Replication/India Networks/MF"+str(village)+".csv"

    leader_filename = data_file+"Matlab Replication/India Networks/HHhasALeader"+str(village)+".csv"
    if os.path.isfile(leader_filename):
        leaders[village] = np.loadtxt(leader_filename, delimiter="\t")

    if os.path.isfile(mf_filename):
        outcome = np.loadtxt(data_file+"Matlab Replication/India Networks/MF"+str(village)+".csv", delimiter=",")
        outcomes[village] = outcome
    print(str(village))

V = len(outcomes.keys())
graphs = {}
datas = [{}]*V
for vil_num, village in enumerate(outcomes.keys()):
    print(round(100*vil_num / V, 1))
    data = copy.deepcopy(As[village])
    graph = nx.to_networkx_graph(data)
    graphs[village] = copy.deepcopy(graph)

    HHIDs_row = {j: i for i, j in enumerate(HHIDs[village])}
    PIDs_row = {j: i for i, j in enumerate(PIDs[village])}

    all_outcomes = [int(outcomes[village][HHIDs_row[i]]) for i in (PIDs[village]-(village*100000))//100]

    sexes = df[df["village"]==village].set_index("pid")["resp_gend"].to_dict()
    ages = df[df["village"]==village].set_index("pid")["age"].to_dict()
    SHGs = df[df["village"]==village].set_index("pid")["shgparticipate"].to_dict()
    SHG_num = {"Yes":1, "No": 0, -999.0: -1}
    sex = [-1 for i in PIDs[village]]
    age = [-1 for i in PIDs[village]]
    SHG = [-1 for i in PIDs[village]]
    for key in sexes.keys():
        sex[PIDs_row[key]] = sexes[key]
        age[PIDs_row[key]] = ages[key]
        SHG[PIDs_row[key]] = SHG_num[SHGs[key]]

    leaders_dict = {i: -1 for i in HHIDs[village]}
    for i in leaders[village]:
        leaders_dict[int(i[0])] = int(i[1])
    is_leader = [leaders_dict[i] for i in (PIDs[village]-(village*100000))//100]

    components = list(len(i) for i in nx.connected_components(graph))

    leads = {i for i, j in enumerate(is_leader) if j == 1}
    mins = {node: 0 for node in graph.nodes()}
    sums = {node: 0 for node in graph.nodes()}
    for node in graph.nodes():
        paths = [nx.shortest_path_length(graph,node,neighbor) for neighbor in set(nx.node_connected_component(graph, node)).intersection(leads).difference(set([node]))]
        if len(paths) == 0:
            mins[node] = 0
            sums[node] = 0
        if len(paths) > 0:
            mins[node] = float(1) / min(paths)
            sums[node] = sum([float(1)/path for path in paths])

    datas[vil_num] = pd.DataFrame({
        "Village":              village,
        "HHID":                 dict(enumerate(PIDs[village]//100)),
        "Outcome":              dict(enumerate(all_outcomes)),
        "Is_Leader":            dict(enumerate(is_leader)),
        "Self_Help":            dict(enumerate(SHG)),
        "Sex":                  dict(enumerate(sex)),
        "Age":                  dict(enumerate(age)),
        "Degree":               graphs[village].degree(),
        "Mean_Neighbor_Degree": nx.average_neighbor_degree(graph),
        "Assortativity":        nx.degree_assortativity_coefficient(graph),
        "LCC_Size":             max(components),
        "Mean_Component_Size":  np.mean(components),
        "Number_Of_Components": len(components),
        "Node_Component_Size":  {node: len(nx.node_connected_component(graph, node)) for node in graph.nodes()},
        "Total_Neighbor_Seeds": {node: len(set(graph.neighbors(node)).intersection(leads)) for node in graph.nodes()},
        "Total_Cluster_Seeds":  {node: len(set(nx.node_connected_component(graph, node)).intersection(leads)) for node in graph.nodes()},
        "Mins":                 mins, # inverse of the minimum path length to an infected node at baseline.
        "Sums":                 sums  # sum of the inverse path length to an infected node at baseline.
    })
    datas[vil_num] = datas[vil_num][[
        "Village","HHID", "Outcome", "Is_Leader", "Self_Help", "Sex", "Age",
        "Degree", "Mean_Neighbor_Degree", "Assortativity",
        "LCC_Size", "Mean_Component_Size", "Number_Of_Components", "Node_Component_Size",
        "Total_Neighbor_Seeds", "Total_Cluster_Seeds", "Mins", "Sums"
    ]]
total_data = pd.concat(datas)
total_data.to_csv("incomplete_data.txt", sep = "\t", index = False)











