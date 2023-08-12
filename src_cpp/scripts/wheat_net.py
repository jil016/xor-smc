import numpy as np
import os

def gen_wheat_data(net_folder):
    # We need:
    # - produces
    # - demand
    # - cost
    # - budget
    # - capacity
    # - disaster models
    if not os.path.exists(net_folder):
        os.mkdir(net_folder)

    wheat = [i for i in range(9)]  # 9
    flour = [i for i in range(9, 16)]  # 7
    bread = [i for i in range(16, 25)]  # 9
    market = [i for i in range(25, 25+19)]  #19
    n_disaster = 4

    # produce
    produce = ['0'] * len(wheat) + ['1'] * len(flour) + ['2'] * len(bread) + ['3'] * len(market)
    with open(net_folder + "/produce.txt", "w") as fp:
        fp.write(f"{market[-1] + 1} {4}\n")
        fp.write(" ".join(produce))

    # demand
    demand = []
    for f in flour:
        demand.append([f, 0, 1])
    for b in flour:
        demand.append([b, 1, 1])
    for m in flour:
        demand.append([m, 2, 1])
    # only need one unit?
    # increase need for 6/7, will increase the production
    with open(net_folder + "/demand.txt", "w") as fp:
        fp.write(f"{len(demand)}\n")
        for line in demand:
            fp.write(" ".join([str(l) for l in line]))
            fp.write("\n")

    # cost
    cost_wf =  [[967,	659,	1330,	810,	1159,	1081,	375],
                [280,	1523,	1637,	1038,	911,	1276,	879],
                [1032,	368,	1042,	1152,	470,	1390,	1212],
                [1623,	1325,	999,	1173,	924,	873,	1543],
                [414,	1112,	1288,	653,	838,	953,	320],
                [1653,	619,	485,	975,	1284,	657,	1327],
                [329,	949,	1125,	464,	649,	736,	263],
                [1143,	120,	680,	485,	950,	350,	830],
                [731,	1181,	1273,	600,	370,	911,	911]]

    cost_fb =  [[1553,	1909,	1353,	900,	757,	450,	1623,	505,	1300],
                [490,	1123,	571,	659,	485,	881,	1325,	773,	100],
                [485,	552,	100,	1230,	661,	1273,	550,	949,	571],
                [975,	1213,	661,	745,	300,	320,	1173,	288,	585],
                [1084,	1378,	1073,	886,	396,	270,	924,	507,	1181],
                [557,	913,	361,	981,	300,	527,	673,	588,	575],
                [1227,	1583,	1031,	300,	370,	711,	1543,	206,	855]]

    cost_bm =  [[743,	866,	720,	520,	200,	485,	300,	400,	300,	619,	320,	400,	1087,	1037,	835,	770,	557,	750,	827],
                [325,	170,	200,	352,	652,	360,	890,	680,	823,	1123,	723,	870,	1736,	1736,	1436,	1335,	913,	1293,	1427],
                [529,	550,	529,	300,	200,	280,	685,	670,	400,	671,	460,	500,	1284,	1184,	1070,	950,	361,	745,	775],
                [1759,	1780,	1765,	1530,	1030,	1530,	978,	1587,	800,	510,	870,	800,	100,	210,	120,	200,	1081,	433,	485],
                [1190,	1363,	1370,	861,	561,	850,	975,	1175,	600,	485,	661,	680,	700,	870,	745,	870,	300,	299,	500],
                [1453,	1673,	1810,	1473,	1273,	1540,	1384,	1680,	1010,	1220,	1160,	710,	760,	904,	886,	986,	727,	810,	976],
                [350,	480,	780,	360,	700,	580,	1313,	1113,	925,	1225,	960,	1100,	1640,	1808,	1618,	1650,	873,	1405,	1200],
                [1478,	1650,	1700,	1100,	949,	1120,	1045,	1525,	905,	773,	995,	960,	300,	701,	581,	650,	588,	587,	868],
                [1100,	1000,	1050,	700,	420,	670,	519,	819,	1855,	180,	250,	200,	730,	670,	620,	559,    425,	200,	304]]

    
    cost_matrix = np.zeros((market[-1] + 1, market[-1] + 1),dtype=int).tolist()
    for i, w_idx in enumerate(wheat):
        for j, f_idx in enumerate(flour):
            cost_matrix[w_idx][f_idx] = cost_wf[i][j]

    for i, f_idx in enumerate(flour):
        for j, b_idx in enumerate(bread):
            cost_matrix[f_idx][b_idx] = cost_fb[i][j]

    for i, b_idx in enumerate(bread):
        for j, m_idx in enumerate(market):
            cost_matrix[b_idx][m_idx] = cost_bm[i][j]

    with open(net_folder + "/cost.txt", "w") as fp:
        for line in cost_matrix:
            fp.write(" ".join([str(l) for l in line]))
            fp.write("\n")
    
    # budget
    # node 0 - 
    budget = ['3000'] * (market[-1] + 1)
    with open(net_folder + "/budget.txt", "w") as fp:
        fp.write(" ".join(budget))
    
    # capacity of each edge -- limited by factory capacity and transportation
    # raw data, then discretized to 0~16
    # should correspond to demand
    # real capacity: 5000, 9000, 8500
    # real demand for market: 2900 bread 
    # real demand for bread factory: 2230 flour
    # real demand for flour factory: 2686 wheat
    # 
    # sum cap * disaster_prob * I >= demand
    # demand / cap * 16^n_disasters
    print("demand = ", np.log2(2900/8500/16 *(16**n_disaster)))

    # rate_wheat_flour = 0.83
    # rate_flour_bread = 1.3
    # rate_market = 1.0
    # probability: discretized by 1/16

    capacity_wf = np.random.randint(10,16, size=(9,7)).tolist()

    capacity_fb = np.random.randint(9,15, size=(7,9)).tolist()

    capacity_bm = np.random.randint(8,14, size=(9,19)).tolist()

    capacity_matrix = np.zeros((market[-1] + 1, market[-1] + 1),dtype=int).tolist()
    for i, w_idx in enumerate(wheat):
        for j, f_idx in enumerate(flour):
            capacity_matrix[w_idx][f_idx] = capacity_wf[i][j]

    for i, f_idx in enumerate(flour):
        for j, b_idx in enumerate(bread):
            capacity_matrix[f_idx][b_idx] = capacity_fb[i][j]

    for i, b_idx in enumerate(bread):
        for j, m_idx in enumerate(market):
            capacity_matrix[b_idx][m_idx] = capacity_bm[i][j]
    
    # print(capacity_matrix)
    # export capacity matrix
    with open(net_folder + "/capacity.txt", "w") as fp:
        fp.write("4\n")
        for line in capacity_matrix:
            fp.write(" ".join([str(l) for l in line]))
            fp.write("\n")

    disaster_models = []
    for i in range(n_disaster):
        disaster_i = np.random.randint(0,2,size=(8,8))
        disaster_models.append(disaster_i.tolist())

    # export disasters
    for i, model in enumerate(disaster_models):
        with open(net_folder + f"/disaster{i}.txt", "w") as fp:
            fp.write("4\n")
            for line in model:
                fp.write(" ".join([str(l) for l in line]))
                fp.write("\n")


    return


def gen_downsized_data(net_folder):
    # We need:
    # - produces
    # - demand
    # - cost
    # - budget
    # - capacity
    # - disaster models
    if not os.path.exists(net_folder):
        os.mkdir(net_folder)

    wheat = [i for i in range(2)]  # 9
    flour = [i for i in range(2, 4)]  # 7
    bread = [i for i in range(4, 6)]  # 9
    market = [i for i in range(6, 8)]  #19
    n_disaster = 4

    # produce
    produce = ['0'] * len(wheat) + ['1'] * len(flour) + ['2'] * len(bread) + ['3'] * len(market)
    
    with open(net_folder + "/produce.txt", "w") as fp:
        fp.write(f"{market[-1] + 1} {4}\n")
        fp.write(" ".join(produce))
    
    # demand
    demand = [[2, 0, 1],
              [3, 0, 1],
              [4, 1, 1],
              [5, 1, 1],
              [6, 2, 1],
              [7, 2, 1]]
    # only need one unit?
    # increase need for 6/7, will increase the production
    with open(net_folder + "/demand.txt", "w") as fp:
        fp.write(f"{len(demand)}\n")
        for line in demand:
            fp.write(" ".join([str(l) for l in line]))
            fp.write("\n")

    # cost
    cost_wf =  [[967,	659],
                [280,	1523]]

    cost_fb =  [[1553,	1909],
                [490,	1123]]
    
    cost_bm =  [[743,	866],
                [325,	170]]
    
    cost_matrix = np.zeros((market[-1] + 1, market[-1] + 1),dtype=int).tolist()
    for i, w_idx in enumerate(wheat):
        for j, f_idx in enumerate(flour):
            cost_matrix[w_idx][f_idx] = cost_wf[i][j]

    for i, f_idx in enumerate(flour):
        for j, b_idx in enumerate(bread):
            cost_matrix[f_idx][b_idx] = cost_fb[i][j]

    for i, b_idx in enumerate(bread):
        for j, m_idx in enumerate(market):
            cost_matrix[b_idx][m_idx] = cost_bm[i][j]

    with open(net_folder + "/cost.txt", "w") as fp:
        for line in cost_matrix:
            fp.write(" ".join([str(l) for l in line]))
            fp.write("\n")
    
    # budget
    # node 0 - 7
    budget = ['3000'] * 8
    with open(net_folder + "/budget.txt", "w") as fp:
        fp.write(" ".join(budget))
    
    # capacity of each edge -- limited by factory capacity and transportation
    # raw data, then discretized to 0~16
    # should correspond to demand
    # real capacity: 5000, 9000, 8500
    # real demand for market: 2900 bread 
    # real demand for bread factory: 2230 flour
    # real demand for flour factory: 2686 wheat
    # 
    # sum cap * disaster_prob * I >= demand
    # demand / cap * 16^n_disasters
    print("demand = ", np.log2(2900/8500/16 *(16**n_disaster)))

    # rate_wheat_flour = 0.83
    # rate_flour_bread = 1.3
    # rate_market = 1.0
    # probability: discretized by 1/16

    capacity_wf = np.random.randint(10,16, size=(2,2)).tolist()

    capacity_fb = np.random.randint(9,15, size=(2,2)).tolist()

    capacity_bm = np.random.randint(8,14, size=(2,2)).tolist()

    capacity_matrix = np.zeros((market[-1] + 1, market[-1] + 1),dtype=int).tolist()
    for i, w_idx in enumerate(wheat):
        for j, f_idx in enumerate(flour):
            capacity_matrix[w_idx][f_idx] = capacity_wf[i][j]

    for i, f_idx in enumerate(flour):
        for j, b_idx in enumerate(bread):
            capacity_matrix[f_idx][b_idx] = capacity_fb[i][j]

    for i, b_idx in enumerate(bread):
        for j, m_idx in enumerate(market):
            capacity_matrix[b_idx][m_idx] = capacity_bm[i][j]
    
    # print(capacity_matrix)
    # export capacity matrix
    with open(net_folder + "/capacity.txt", "w") as fp:
        fp.write("4\n")
        for line in capacity_matrix:
            fp.write(" ".join([str(l) for l in line]))
            fp.write("\n")

    disaster_models = []
    for i in range(n_disaster):
        disaster_i = np.random.randint(0,2,size=(8,8))
        disaster_models.append(disaster_i.tolist())

    # export disasters
    for i, model in enumerate(disaster_models):
        with open(net_folder + f"/disaster{i}.txt", "w") as fp:
            fp.write("4\n")
            for line in model:
                fp.write(" ".join([str(l) for l in line]))
                fp.write("\n")

    return

if __name__ == '__main__':
    np.random.seed(1086)
    gen_wheat_data("./test_net_large")
    gen_downsized_data("./test_net")