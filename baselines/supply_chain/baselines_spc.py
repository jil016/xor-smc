import numpy as np
from scipy.optimize import linprog
# MCMC greedy
from supply_chain_loader import SupplyNet


class task(object):
    def __init__(self, filename):
        self.supply_chain = SupplyNet(filename)

    def ground_truth_production(self, ):
        res = {}
        cost_dict = {}
        for my_node in range(self.N):
            # inspect every in-edge
            expected_capacities = [0] * self.N
            my_budget = self.budget[my_node]
            my_cost = [0] * self.N

            for s in range(self.N):
                if self.capacity[s][my_node] != 0:

                    # estimate the expectation
                    # ...
                    # linear programming to get a selection
                    prob_list = []
                    for dis in self.disasters:
                        prob_list.append(dis[s][my_node] / 16.0)

                    prob_pass = 1.0
                    for p in prob_list:
                        prob_pass = prob_pass * (1 - p)
                    # calc expected utility
                    expected_capacities[s] = prob_pass * self.capacity[s][my_node]
                    my_cost[s] = self.cost[s][my_node]
            res[my_node] = expected_capacities
            cost_dict[my_node] = my_cost
            # print(my_cost)
            # print(expected_capacities)
            # consider my_budget

        # calc GT best expected production using LP
        return res, cost_dict


def linearProgramming(exp_cap, cost, budget=3000):
    c = [-i for i in exp_cap]
    A = [cost]
    b = [budget]
    x_bounds = [(0, 1)] * len(A)
    res = linprog(c, A_ub=A, b_ub=b, bounds=x_bounds)
    print(res.fun)
    print(res.x)
    return res.fun


def mixedIntegerProgramming(exp_cap, cost, budget=3000):
    from mip import Model, xsum, maximize, BINARY

    p = exp_cap
    w = cost
    c, I = budget, range(len(w))

    m = Model()

    x = [m.add_var(var_type=BINARY) for i in I]
    m.objective = maximize(xsum(p[i] * x[i] for i in I))
    m += xsum(w[i] * x[i] for i in I) <= c
    m.optimize()

    selected = [i for i in I if x[i].x >= 0.99]
    print("selected items: {}".format(selected))
    return


def mixedIntegerProgrammingSPC(model, trades, budget):
    return


def evaluateXORres():
    N = 100
    file_path = "./LOG-SPC-test/result.log"
    network_folder = "./test_net"

    # read XOR result
    res = []
    with open(file_path, "r") as fp:
        fp.readline()
        for i in range(N):
            line = fp.readline()
            line = line[:-2].split(' ')
            row = [int(l) for l in line]
            res.append(row)

    # read GT
    supply_net = SupplyNet(network_folder, 10)
    res_dict, cost_dict = supply_net.ground_truth_production()

    # calc XOR quality
    xor_res = []
    for e in supply_net.end_users:
        cap_acc = 0
        for i in range(N):
            if res[i][e] == 1:
                cap_acc += res_dict[e][i]
        xor_res.append(cap_acc)

    print(supply_net.end_users)
    print(xor_res)

    gt_res = []
    for e in supply_net.end_users:
        # temp = mixedIntegerProgramming(res_dict[e], cost_dict[e])
        # gt_res.append(-temp)
        best = 0
        # simply pick best
        for i, c in enumerate(cost_dict[e]):
            if c <= 1000:
                if res_dict[e][i] > best:
                    best = res_dict[e][i]

        gt_res.append(best)
        # print(">>>>>>>>>>>>>>>>>>>>")
        # print(res_dict[e])
        # print("====================")
        # print(cost_dict[e])

    print(gt_res)
    return


if __name__ == '__main__':
    evaluateXORres()
    # sample a trading plan from the UAI file
