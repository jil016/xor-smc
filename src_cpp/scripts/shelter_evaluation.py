import numpy as np
import os

def generateCNF(outfolder, sources, sinks):
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)

    
    return




if __name__ == '__main__':

    sources = [0, 10, 20]

    xor_121 = [19, 20, 43, 77, 105]
    xor_183 = [0, 20, 24, 121, 157]
    xor_246 = [0, 10, 18, 119, 147]
    xor_388 = [1, 13, 39, 87, 175]

    gibbs_121 = []
    gibbs_183 = []
    gibbs_246 = []
    gibbs_388 = []

    unigen_121 = []
    unigen_183 = []
    unigen_246 = []
    unigen_388 = []

    quick_121 = []
    quick_183 = []
    quick_246 = []
    quick_388 = []

    generateCNF("shelter_CNFs_121", sources, sinks)