# Solving Satisfiability Modulo Counting for Symbolic and Statistical AI Integration With Provable Guarantees






# Generate the data for supply chain




# Install Dependency
- python related packages:
```cmd
pip install -r requirements.txt
```
- C++ related packages.  Our algorithm use the [CPLEX solver](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer), so we need to install the CPLEX before running the algorithm.


# Useful commands for CPLEX

```
sudo chown -R username /ibm_directory/
```

# Applications
Please goto the `src` folder to run the program for our XOR-SMC the two applications.
- shelter design.
- supply chain.

the baselines algorithms are collected in the folder with name `baselines`.

# Cite

```
@article{DBLP:proceedings/AAAI24/li-xor-smc,
  author       = {Jinzhao Li and
                  Nan Jiang and
                  Yexiang Xue},
  title        = {Solving Satisfiability Modulo Counting for Symbolic and Statistical
                  {AI} Integration With Provable Guarantees},
  journal      = {AAAI},
  year         = {2024}
}
```