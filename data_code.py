# ############################################################################################ #
# #############                                                           #################### #
# #############       Infectious Simulation in Arbitrary Networks         #################### #
# #############                                                           #################### #
# ############################################################################################ #

import os, sys, random, math, copy
import numpy as np, networkx as nx#, lifelines
from string import *
from operator import itemgetter
import pandas as pd
#exec(open("C:\\Phoenix\\School\\Miscellany\\core.py").read())
#os.chdir("C:\\Users\\Patrick\\Desktop")

trial_num = sys.argv[1]
setup = sys.argv[2]

if(setup[0] == "0"):
    setup_high_degree = "No"

if(setup[0] == "1"):
    setup_high_degree = "Yes"

if(setup[1] == "0"):
    setup_degree_distribution = "Poisson"

if(setup[1] == "1"):
    setup_degree_distribution = "Powerlaw"

if(setup[2] == "0"):
    setup_assortativity = "Negative"

if(setup[2] == "1"):
    setup_assortativity = "Positive"

if(setup[3] == "0"):
    setup_block_structure = "No"

if(setup[3] == "1"):
    setup_block_structure = "Yes"

if(setup[4] == "0"):
    setup_infectivity = "Unit"

if(setup[4] == "1"):
    setup_infectivity = "Degree"

if(setup[5] == "0"):
    setup_baseline_prevalence= "Low"

if(setup[5] == "1"):
    setup_baseline_prevalence= "High"



data_filepath = "Data"+setup
if not os.path.exists("Data"+setup):
    os.makedirs("Data"+setup)

np.set_printoptions(precision=3, suppress = True)

def logit(x):
    return np.log(x) - np.log(1-x)

def expit(x):
    return float(1) / (1 + math.exp(-x))

def prob_func(x):
    return expit(x*1/2) # after standardizing the obs variable, should give close to 10% flip at -3, and 90% flip at 3.

def weighted_sample(choices, n):
    total   = sum(weight for value, weight in choices)
    randoms = sorted([random.uniform(0, total) for r in range(n)])
    choices = sorted(choices, key=itemgetter(1),reverse=True)
    sample  = []
    weight_sum = 0
    r = 0
    for choice, weight in choices:
        weight_sum += weight
        if r >= n:
            break
        if weight_sum > randoms[r]:
            sample.append(choice)
            r += 1
    return sample

def find_assortativity(graph):
    j = np.array([graph.degree(i[0]) for i in graph.edges()]) - 1
    k = np.array([graph.degree(i[1]) for i in graph.edges()]) - 1
    M = len(graph.edges())
    Minv = float(1) / M
    top = (Minv*sum(j*k)-(Minv/2*sum(j+k))**2)
    bottom = (Minv/2*sum(j**2+k**2)-(Minv/2*sum(j+k))**2)
    return(float(top) / bottom)

obs_dependence_strength = 1 # new
trials = 1
n = 120
runs = 2#1000
clusters = 102
ns = [120]*(clusters//6)+[200]*(clusters//6)+[280]*(clusters//6)+[120]*(clusters//6)+[200]*(clusters//6)+[280]*(clusters//6)#[200]*clusters
assortative_cluster = [True]*(clusters//2) + [True]*(clusters//2)#probably delete

overall_average_degree = 1
if(setup_high_degree == "Yes"):
    overall_average_degree = 5

average_degree = [overall_average_degree]*(clusters//2) + [overall_average_degree]*(clusters//2)
C = 4

if setup_block_structure == "No":
    Lambda = .0
    Mu = .0

if setup_block_structure == "Yes":
    Lambda = .1
    Mu = .8

if setup_assortativity == "Negative":
    assortativity_threshold = -0.3 # +-.3, .01

if setup_assortativity == "Positive":
    assortativity_threshold = 0.3 # +-.3, .01

M = 5  # number of epidemics to run.  Just in case, this is more than one.
#prop_seeds = 1#5
ps = [.30, .10]
IV_ps = np.array([.30, .25])*0 # IV quantile is inherent here.  range from 0 to 1.
quantile = 25 #+ random.choice(range(-4,4))

if setup_baseline_prevalence == "Low":
    prop_seeds = 1
    baseline_quantile  = 2

if setup_baseline_prevalence == "High":
    prop_seeds = 10
    baseline_quantile  = 25

study_end_incidence = 5
study_end = 5 # made to represent years

baseline_total = float(baseline_quantile)/100 * sum(ns)
study_end_total = baseline_total + float(study_end_incidence)/100 * (1-float(baseline_quantile)/100) * sum(ns)  # this adds the incidence of the number of individuals available after baseline.

if setup_infectivity == "Unit":
    infectivity = ["unit", "unit"]
    concurrency = ["unit","unit"]

if setup_infectivity == "Degree":
    infectivity = ["degree", "degree"]
    concurrency = ["degree","degree"]

verbose = False
timeout = 20

for trial in [trial_num]:#range(trials)
    total_data = pd.DataFrame({})
    graphs = [nx.Graph()]*clusters
    gs = [{}]*clusters
    adjusted_ps = [{}]*clusters
    ccs = [{}]*clusters # concurrency per node
    iis = [{}]*clusters # infectivity per node
    LCCs = [{}]*clusters
    seedss = [{}]*clusters
    degreess = [{}]*clusters
    infecteds = [{}]*clusters
    ts = [0]*M

    # Generate a bipartite degree-corrected stochastic blockmodel with assortative rewiring (preserving degree and block)
    for cluster in range(clusters):
        while len(LCCs[cluster]) < float(baseline_quantile+study_end_incidence)/100* ns[cluster]:  # this ensures that every graph has a LCC that can sustain a sufficiently-sized epidemic.
            assortative = True#assortative_cluster[cluster]
            if setup_degree_distribution == "Poisson":
                k_female = np.random.poisson(average_degree[cluster],ns[cluster]/2)
                k_male   = np.random.poisson(average_degree[cluster],ns[cluster]/2)
            if setup_degree_distribution == "Powerlaw":
                k_female =  ((average_degree[cluster]-1)*np.random.uniform(size=ns[cluster]/2)+1).astype(int)*np.random.zipf(2.5, size=ns[cluster]/2) # will require some editing if mean(degree) differs from K.
                k_male   =  ((average_degree[cluster]-1)*np.random.uniform(size=ns[cluster]/2)+1).astype(int)*np.random.zipf(2.5, size=ns[cluster]/2)
            
            counter_threshold = 30
            
            k = np.concatenate((k_female, k_male)) # probably just have to have the same sum (maybe not even that), but this assumes they have the same distribution.  For now!
            gs[cluster] = {i: 2*C*i//ns[cluster] for i in range(ns[cluster])}
            kappa = [np.sum([k[i] for i in range(ns[cluster]) if gs[cluster][i] == K]) for K in range(2*C)]
            m = float(sum(kappa))/2
            theta = [float(k[i]) / kappa[gs[cluster][i]] for i in range(ns[cluster])]
            omega_random = np.zeros((2*C,2*C))
            omega_zeros = np.zeros((C,C))
            omega_block = np.zeros((C,C))
            for i in range(C):
                for j in range(C):
                    omega_block[i,j] = float(kappa[i]*kappa[j]) / (2*m)
                omega_random[:C,C:(2*C)] = omega_block
                omega_random[C:(2*C),:C] = omega_block
            
            omega_planted = np.zeros((2*C,2*C))
            for i in range(2*C):
                for j in range(2*C):
                    if i == ((j+C) % (2*C)):
                        omega_planted[i,j] = float(kappa[i])/2
            
            omega_sexwork = copy.deepcopy(omega_planted)
            omega_sexwork[:C,C] = copy.deepcopy(omega_random[:C,C])
            omega_sexwork[C,:C] = copy.deepcopy(omega_random[C,:C])
            omega_sexwork[C:,0] = copy.deepcopy(omega_random[C:,0])
            omega_sexwork[0,C:] = copy.deepcopy(omega_random[0,C:])
            
            omega = Lambda * omega_planted + Mu * omega_sexwork + (1 - Lambda - Mu) * omega_random
            
            graph = nx.Graph()
            graph.add_nodes_from(range(ns[cluster]))
            for i in range(2*C):
                for j in range(2*C):
                    mean_edges = int(omega[i,j])
                    if i == j: #True
                        mean_edges = float(mean_edges) / 2
                    pois = np.random.poisson(mean_edges)
                    ns_i = [node for node in range(ns[cluster]) if gs[cluster][node] == i]
                    ns_j = [node for node in range(ns[cluster]) if gs[cluster][node] == j]
                    p_i = [theta[node] for node in ns_i]
                    p_j = [theta[node] for node in ns_j]
                    stubs_i = [np.random.choice(ns_i,p=p_i) for stub in range(pois)]
                    stubs_j = [np.random.choice(ns_j,p=p_j) for stub in range(pois)]
                    graph.add_edges_from(zip(stubs_i, stubs_j))
            LCCs[cluster] = sorted(nx.connected_components(graph), key = len, reverse = True)[0]
            
            degrees = list(graph.degree().values())
            m = len(graph.edges())
            
            counter       = 0
            assortativity = 0
            edges = graph.edges()
            omega_edges = [[[edge for edge in edges if ((gs[cluster][edge[0]] == i and gs[cluster][edge[1]] == j) or (gs[cluster][edge[0]] == j and gs[cluster][edge[1]] == i))] for j in range(2*C)] for i in range(2*C)] # a matrix of edge lists according to 2*C membership.  Usage: omega_edges[i][j] is a bunch of edges.
            #
            while assortativity_threshold - ((2*int(assortative)-1) * assortativity) > 0 and counter < counter_threshold:
                for i in range(100):
                    random_edge = random.choice(graph.edges())
                    same_block_edge = random.choice(omega_edges[gs[cluster][random_edge[0]]][gs[cluster][random_edge[1]]]) # select another edge from the same omega block.
                    current_edges = (random_edge, same_block_edge)
                    remaining_degrees = [[degrees[stub] for stub in edge] for edge in current_edges]
                    current_degree_difference  = abs(remaining_degrees[0][0] - remaining_degrees[0][1]) + abs(remaining_degrees[1][0] - remaining_degrees[1][1])
                    rewiring_degree_difference = abs(remaining_degrees[0][0] - remaining_degrees[1][1]) + abs(remaining_degrees[0][1] - remaining_degrees[1][0])
                    if (2*int(assortative)-1) * (current_degree_difference - rewiring_degree_difference) > 0 :
                        new_edges = ((random_edge[0],same_block_edge[1]),(random_edge[1],same_block_edge[0]))
                        graph.remove_edges_from(current_edges)
                        graph.add_edges_from(new_edges)
                assortativity = find_assortativity(graph)
                counter += 1
        graphs[cluster] = graph
    
    # end network generation.  start the epidemic.
    for cluster in range(clusters):
        seedss[cluster] = random.sample(graphs[cluster].nodes(), int(max(float(n)*float(prop_seeds)/100,1))) #LCCs[cluster]
        degreess[cluster] = list(graphs[cluster].degree().values())
        infecteds[cluster] = set(seedss[cluster])
        ccs[cluster] = {node: (graphs[cluster].degree(node) - 1)*(int(concurrency[node//ns[cluster]] == "degree")) + 1 for node in range(ns[cluster])}
        iis[cluster] = {node: (graphs[cluster].degree(node) - 1)*(int(infectivity[node//ns[cluster]] == "degree")) + 1 for node in range(ns[cluster])}
        adjusted_ps[cluster] = {node: ps[0]*float(iis[cluster][node])/(ccs[cluster][node] + int(ccs[cluster][node]==0)) for node in range(ns[cluster])}
        t = 0
        if False:   #cluster == 0 is better. this probably has to be restructured.
            print("trial " + str(trial))
            print("proportion of target degree:" + str(float(sum(degrees)) / sum(k))) #always slightly low...
            print("assort: " + str(find_assortativity(graph)))
            print(assortativity, counter)
            print(np.mean(list(graph.degree().values())), np.mean(k))
            print(float(len(seedss[cluster]))/ns[cluster])
    
    current_prevalence = sum(len(infecteds[cluster]) for cluster in range(clusters))
    class Baseline_Reached(BaseException): pass
    try:
        for t in range(timeout):
            for cluster in random.sample(range(clusters), clusters):
                newly_infected = set()
                for infected_node in random.sample(infecteds[cluster],len(infecteds[cluster])):
                    neighbors = graphs[cluster].neighbors(infected_node)
                    if len(neighbors) > 0:
                        for neighbor in random.sample(neighbors, ccs[cluster][infected_node]):
                            successful_infection = random.random() < adjusted_ps[cluster][infected_node]
                            if successful_infection and neighbor not in infecteds[cluster]:
                                current_prevalence += 1
                                newly_infected.add(neighbor)
                                infecteds[cluster] = infecteds[cluster].union(newly_infected)
                                if current_prevalence >= baseline_total:
                                    raise Baseline_Reached()
    except Baseline_Reached:
        pass

    prevalences = copy.deepcopy(infecteds)
    prevalence_outcomes = [{}]*clusters
    for cluster in range(clusters):
        prevalence_outcomes[cluster] = copy.deepcopy({node: 0 for node in graphs[cluster].nodes()})
        for infected_node in prevalences[cluster]:
            prevalence_outcomes[cluster][infected_node] = 1

    obs_file = data_filepath + "/" + "obs_parameters.txt"
    sum_neighbor_seeds = {cluster: sum([len(set(graphs[cluster].neighbors(node)).intersection(set(prevalences[cluster]))) for node in range(ns[cluster])]) for cluster in range(clusters)}
    obs = copy.deepcopy(sum_neighbor_seeds)
    if not os.path.isfile(obs_file):
        pd.DataFrame([np.mean(list(obs.values())), np.std(list(obs.values()))]).to_csv(obs_file,index=False, header=False)
    params = np.loadtxt(obs_file)
    obs = {cluster: (obs[cluster]-params[0])/params[1] for cluster in range(clusters)}
    flip_probs = {cluster: prob_func(obs[cluster]*obs_dependence_strength) for cluster in range(clusters)}
    trts = {cluster: np.random.binomial(1, flip_probs[cluster]) for cluster in range(clusters)}



    counter_inf = np.zeros((runs*2,clusters))
    for i in range(runs*2+1):
        if i < runs:
            counter_ps = (ps[0],ps[0])
        if i > runs:
            counter_ps = (ps[1],ps[1])
        if i == 2*runs:
            counter_ps = copy.deepcopy(ps)
        infecteds = copy.deepcopy(prevalences)
        for cluster in range(clusters):
            current_incidence = sum(len(infecteds[cluster]) for cluster in range(clusters))
            adjusted_ps[cluster] = {node: counter_ps[trts[cluster]]*float(iis[cluster][node])/(ccs[cluster][node] + int(ccs[cluster][node]==0)) for node in range(ns[cluster])}

        class StudyEnd_Reached(BaseException): pass
        try:
            for t in range(study_end):#range(timeout):
                for cluster in random.sample(range(clusters), clusters):
                    newly_infected = set()
                    for infected_node in random.sample(infecteds[cluster],len(infecteds[cluster])):
                        neighbors = graphs[cluster].neighbors(infected_node)
                        if len(neighbors) > 0:
                            for neighbor in random.sample(neighbors, ccs[cluster][infected_node]):
                                successful_infection = random.random() < adjusted_ps[cluster][infected_node]
                                if successful_infection and neighbor not in infecteds[cluster]:
                                    current_incidence += 1
                                    newly_infected.add(neighbor)
                                    infecteds[cluster] = infecteds[cluster].union(newly_infected)
                                    if current_incidence >= study_end_total: # t>study_end:
                                        raise StudyEnd_Reached()
        except StudyEnd_Reached:
            pass
        if i < 2*runs:
            for j in range(clusters):
                counter_inf[i][j] = len(infecteds[j].difference(prevalences[j]))

    means_0 = np.mean(counter_inf[range(runs),       :]/np.array(ns),axis=0)
    means_A = np.mean(counter_inf[range(runs,2*runs),:]/np.array(ns),axis=0)
    cluster_prevalences = np.array([len(i) for i in prevalences])/np.array(ns)

    with open(data_filepath + "/" + "logOR_"+str(trial_num)+".txt", "w") as f:
        f.write(
            str(np.mean(means_A+.000001))+"\t"+str(np.mean(means_0+.000001))+"\n"+
            str(np.mean((means_A+.000001)/(1-cluster_prevalences)))+"\t"+str(np.mean((means_0+.000001)/(1-cluster_prevalences)))+"\n"+
            str(np.mean(means_A+cluster_prevalences))+"\t"+str(np.mean(means_0+cluster_prevalences))+"\n"
        )

    #
    datas = [{}]*clusters
    for cluster in range(clusters):
        outcomes = {node: 0 for node in graphs[cluster].nodes()}
        for infected_node in infecteds[cluster]:
            outcomes[infected_node] = 1
        #
        mins = {node: 0 for node in range(ns[cluster])}
        sums = {node: 0 for node in range(ns[cluster])}
        for node in range(ns[cluster]):
            paths = [nx.shortest_path_length(graphs[cluster],node,neighbor) for neighbor in set(nx.node_connected_component(graphs[cluster], node)).intersection(prevalences[cluster]).difference(set([node]))]
            if len(paths) == 0:
                mins[node] = 0
                sums[node] = 0
            if len(paths) > 0:
                mins[node] = float(1) / min(paths)
                sums[node] = sum([float(1)/path for path in paths])
        
        components = [len(C) for C in nx.connected_components(graphs[cluster])]
        datas[cluster] = pd.DataFrame({
            # main informations relevant to basic analysis
            "Trt":                  trts[cluster],
            "Cluster":              cluster,
            "Outcome":              outcomes,#dict(zip(range(n),np.array(list(outcomes.values()))-np.array(list(prevalence_outcomes[cluster].values())))),
            
            # Degree-based covariates
            "Degree":               graphs[cluster].degree(),
            "Mean_Neighbor_Degree": nx.average_neighbor_degree(graphs[cluster]),
            "Assortativity":        nx.degree_assortativity_coefficient(graphs[cluster]),
            
            # Community-based covariates
            "Sex_Worker":           {node: int((node < int(float(ns[cluster])/(C*2))) or (node>(ns[cluster]//2) and node < (ns[cluster]//2 +ns[cluster]//(C*2)))) for node in range(ns[cluster])},
            
            # Component-based covariates
            "LCC_Size":             max(components),
            "Mean_Component_Size":  np.mean(components),
            "Number_Of_Components": len(components),
            "Node_Component_Size":  {node: len(nx.node_connected_component(graphs[cluster], node)) for node in range(ns[cluster])},
            
            # Infectious path-based covariates
            "Total_Neighbor_Seeds": {node: len(set(graphs[cluster].neighbors(node)).intersection(set(prevalences[cluster]))) for node in range(ns[cluster])},
            "Total_Cluster_Seeds":  {node: len(set(nx.node_connected_component(graphs[cluster], node)).intersection(prevalences[cluster])) for node in range(ns[cluster])},
            "Mins":                 mins, # inverse of the minimum path length to an infected node at baseline.
            "Sums":                 sums, # sum of the inverse path length to an infected node at baseline.
            
            # Ancillary network featres not yet being utilized
            "Age_Group":            {node: gs[cluster][node]%C for node in range(ns[cluster])},
            "Sex":                  {node: gs[cluster][node]//C for node in range(ns[cluster])},
            "IV_weights":           0,
            "Prevalences":          prevalence_outcomes[cluster],
            "Obs": obs[cluster]
        })
        
        ##     outcomes[m]                = outcome
        ##      total_neighbor_outcomes[m] = total_neighbor_outcome
        ##      data["total_neighbor_outcomes"+"_"+hypothesis+"_"+str(m)] = pd.Series(total_neighbor_outcome)
        ##      data["Outcome"+"_"+hypothesis+"_"+str(m)] = pd.Series(outcome)
        #
        datas[cluster] = datas[cluster][[
        "Trt", "Cluster", "Outcome",
        "Degree", "Mean_Neighbor_Degree", "Assortativity",
        "Sex_Worker",
        "LCC_Size", "Mean_Component_Size", "Number_Of_Components", "Node_Component_Size",
        "Total_Neighbor_Seeds", "Total_Cluster_Seeds", "Mins", "Sums",
        "Age_Group", "Sex", "IV_weights", "Prevalences", "Obs"
        ]]
    total_data = pd.concat(datas)
    total_data.to_csv(data_filepath + "/" + "CRT_"+str(trial_num)+".txt", sep = "\t", index = False)

print(omega)
for i in omega_edges:
    print([len(j) for j in i])


min(list(flip_probs.values()))
max(list(flip_probs.values()))









