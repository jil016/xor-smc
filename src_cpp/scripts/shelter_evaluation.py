import numpy as np
import os
from utils_data import *
from graph import Graph

def generateCNF(outfolder, sources, sinks, graph_file):
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)

    graph = Graph()
    graph.readFromFile(graph_file)
    N = graph.N
    ctx = bx.Context()
    for src in sources:
        for sink in sinks:
            flow = [[ctx.get_var(f'x_{i}_{j}') for j in range(N)]
                                for i in range(N)]

            cnf, sub_cnfs = flow2CNF_nbf(graph, flow, src, sink)
            exportCNF(cnf, outfolder + f"/src{src}_sink{sink}")
            print(f"/src{src}_sink{sink} exported!")

            # cnt = 0
            # for sat in cnf.iter_sat():
            #     # print(sat)
            #     cnt += 1
            # print(cnt)
    return


def calcDirectly(outfolder, sources, sinks, graph_file):
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)

    graph = Graph()
    graph.readFromFile(graph_file)
    N = graph.N
    ctx = bx.Context()
    for src in sources:
        for sink in sinks:
            flow = [[ctx.get_var(f'x_{i}_{j}') for j in range(N)]
                                for i in range(N)]

            cnf, _ = flow2CNF_nbf(graph, flow, src, sink)
            
            cnt = 0
            for sat in cnf.iter_sat():
                # print(sat)
                cnt += 1
            
            with open(outfolder + f"/src{src}_sink{sink}.res", "w") as fp:
                fp.write(str(cnt))

            print(f"Result written to src{src}_sink{sink}.res!")
    return


def useDirectResults(eval_folder, xor_res, gibbs_res, quick_res, unigen_res):
    sources = [0, 10, 20]

    xor_stats = []
    gibbs_stats = []
    quick_stats = []
    unigen_stats = []

    for src in sources:
        path_cnt = 0
        for sink in xor_res:
            with open(eval_folder + f"/src{src}_sink{sink}.res", "r") as fp:
                num = fp.readline()
                path_cnt += int(num)

        xor_stats.append(path_cnt)

    print(xor_stats)

    for src in sources:
        path_cnt = 0
        for sink in gibbs_res[:5]:
            with open(eval_folder + f"/src{src}_sink{sink}.res", "r") as fp:
                num = fp.readline()
                path_cnt += int(num)

        gibbs_stats.append(path_cnt)

    print(gibbs_stats)


    for src in sources:
        path_cnt = 0
        for sink in quick_res[:5]:
            with open(eval_folder + f"/src{src}_sink{sink}.res", "r") as fp:
                num = fp.readline()
                path_cnt += int(num)

        quick_stats.append(path_cnt)

    print(quick_stats)


    for src in sources:
        path_cnt = 0
        for sink in unigen_res[:5]:
            with open(eval_folder + f"/src{src}_sink{sink}.res", "r") as fp:
                num = fp.readline()
                path_cnt += int(num)

        unigen_stats.append(path_cnt)

    print(unigen_stats)

    return

def useSharpSatResults(eval_folder, xor_list, gibbs_list, quick_list, unigen_list):
    sources = [0, 10, 20]

    xor_stats = []
    gibbs_stats = []
    quick_stats = []
    unigen_stats = []

    for src in sources:
        path_cnt = 0
        for sink in xor_list:
            with open(eval_folder + f"/src{src}_sink{sink}.out", "r") as fp:
                lines = fp.readlines()
                log_cnt = lines[-2][:-1].split(" ")[-1]
                int_cnt = lines[-1][:-1].split(" ")[-1]
                path_cnt += int(int_cnt)
        xor_stats.append(path_cnt)
    print(xor_stats)

    for src in sources:
        path_cnt = 0
        for sink in gibbs_list:
            with open(eval_folder + f"/src{src}_sink{sink}.out", "r") as fp:
                lines = fp.readlines()
                log_cnt = lines[-2][:-1].split(" ")[-1]
                int_cnt = lines[-1][:-1].split(" ")[-1]
                path_cnt += int(int_cnt)
        gibbs_stats.append(path_cnt)
    print(gibbs_stats)

    for src in sources:
        path_cnt = 0
        for sink in quick_list:
            with open(eval_folder + f"/src{src}_sink{sink}.out", "r") as fp:
                lines = fp.readlines()
                log_cnt = lines[-2][:-1].split(" ")[-1]
                int_cnt = lines[-1][:-1].split(" ")[-1]
                path_cnt += int(int_cnt)
        quick_stats.append(path_cnt)
    print(quick_stats)


    for src in sources:
        path_cnt = 0
        for sink in unigen_list:
            with open(eval_folder + f"/src{src}_sink{sink}.out", "r") as fp:
                lines = fp.readlines()
                log_cnt = lines[-2][:-1].split(" ")[-1]
                int_cnt = lines[-1][:-1].split(" ")[-1]
                path_cnt += int(int_cnt)
        unigen_stats.append(path_cnt)
    print(unigen_stats)

    return

# 121
# [972, 1080, 1512]
# [288, 648, 2376]
# [1044, 864, 1620]
# [828, 648, 1512]




# 183
# [38880, 25920, 46656]
# [18144, 10368, 18144]
# [10368, 10368, 25920]
# [10368, 5184, 25920]

# 246
# [229864635, 255537828, 137151144]
# [165337200, 212182740, 151047342]
# [168919506, 219347352, 154117890]
# [168919506, 219347352, 154117890]


# 388
# [770943744, 2448880128, 256981248]
# [1027924992, 2403530496, 64245312]
# [342641664, 1133740800, 64245312]
# [342641664, 1133740800, 64245312]

#hawaii121-xor-gibbs-quick-unigen
#972, 288, 864, 648

#hawaii183-xor-gibbs-quick-unigen
#25920, 10368, 10368, 5184

#hawaii246-xor-gibbs-quick-unigen
#137151144, 151047342, 154117890, 154117890

#hawaii388-xor-gibbs-quick-unigen
#770943744, 64245312, 64245312, 64245312

def mainFunc():
    sources = [0, 10, 20]
    graph121 = "graphs/graph_hawaii_121.txt"
    graph183 = "graphs/graph_hawaii_200.txt"
    graph246 = "graphs/graph_hawaii_250.txt"
    graph388 = "graphs/graph_hawaii_388.txt"

    xor_121 = [19, 20, 43, 77, 105]
    xor_183 = [0, 20, 24, 121, 157]
    xor_246 = [0, 10, 18, 119, 147]
    xor_388 = [1, 13, 39, 87, 175]

    gibbs_121 = [37, 91, 112, 77, 65, 30, 112, 51, 43, 58]
    gibbs_183 = [165, 10, 110, 162, 7, 76, 24, 182, 29, 100]
    gibbs_246 = [166, 10, 218, 242, 172, 32, 200, 113, 222, 240]
    gibbs_388 = [165, 10, 314, 84, 42, 176, 122, 85, 371, 311]

    quick_121 = [37, 41, 90, 114, 45, 120, 32, 72, 113, 94]
    quick_183 = [165, 41, 173, 32, 113, 165, 50, 106, 142, 48]
    quick_246 = [16, 121, 218, 242, 234, 32, 200, 114, 222, 60]
    quick_388 = [165, 10, 218, 173, 248, 32, 200, 113, 165, 50]

    unigen_121 = [37, 107, 90, 17, 45, 120, 32, 72, 113, 94]
    unigen_183 = [165, 10, 173, 32, 113, 165, 50, 106, 142, 48]
    unigen_246 = [11, 10, 219, 222, 172, 32, 102, 203, 222, 239]
    unigen_388 = [165, 10, 218, 173, 248, 32, 200, 113, 165, 50]

    hawaii121_union = xor_121 + gibbs_121 + quick_121 + unigen_121
    hawaii121_union = np.unique(hawaii121_union).tolist()
    # print(hawaii121_union)

    hawaii183_union = xor_183 + gibbs_183 + quick_183 + unigen_183
    hawaii183_union = np.unique(hawaii183_union).tolist()
    # print(hawaii183_union)

    hawaii246_union = xor_246 + gibbs_246 + quick_246 + unigen_246
    hawaii246_union = np.unique(hawaii246_union).tolist()
    print(hawaii246_union)

    hawaii388_union = xor_388 + gibbs_388 + quick_388 + unigen_388
    hawaii388_union = np.unique(hawaii388_union).tolist()
    print(hawaii388_union)

    # calcDirectly("shelter_eval_121", sources, hawaii121_union, graph121)
    # calcDirectly("shelter_eval_183", sources, hawaii183_union, graph183)

    # useDirectResults("shelter_eval_121", xor_121, gibbs_121, quick_121, unigen_121)
    # useDirectResults("shelter_eval_183", xor_183, gibbs_183, quick_183, unigen_183)

    # generateCNF("shelter_CNFs_250", sources, hawaii246_union, graph246)
    # generateCNF("shelter_CNFs_388", sources, hawaii388_union, graph388)
    # for 250 and 388, must use sharpSAT-td as in shelter_eval.sh
    # 
    # Then
    useSharpSatResults("LOG-Eval/shelter_CNFs_250", xor_246, gibbs_246[:5], quick_246[:5], unigen_246[:5])

    useSharpSatResults("LOG-Eval/shelter_CNFs_388", xor_388, gibbs_388[:5], quick_388[:5], unigen_388[:5])

    return

if __name__ == '__main__':
    mainFunc()
