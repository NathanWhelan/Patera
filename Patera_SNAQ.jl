using Distributed
addprocs(10)

using Pkg
using PhyloNetworks
using PhyloPlots
using CSV
using RCall
using QuartetNetworkGoodnessFit, DataFrames, CSV
using Gadfly
cd("C:\\Users\\nwhelan\\Documents\\Patera\\RADseq\\SNAQ\\m3M5n5_update")



iqtreeTrees_BS10=joinpath("patera_m3M5n5_populations_r100_maf025_multiSNP_phylo.ALL-pars.BS10.trees")
genetrees_BS10=readMultiTopology(iqtreeTrees_BS10)
#genetrees[3]  #Sanity check
#plot(genetrees[3], :R) #Sainty Check.
q,t = countquartetsintrees(genetrees_BS10) #Generates quartet values (not needed if done with BUCKy)
df = writeTableCF(q,t) #Converts quartets to table.
CSV.write("tableCF_BS10.csv", df) #writes table for use later
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv") #probably could use df, but wrote to disk and now read back in.
#less("tableCF.csv") #Can be done for sanity check.


astralfile_rooted = ("patera_m3M5n5_populations_r100_maf025_multiSNP_phylo.ALL-pars.astral.rooted.tre")
astraltree_rooted = readTopology(astralfile_rooted)
#Read in Topology Level 1 for starting Tree
T_rooted=readTopologyLevel1("patera_m3M5n5_populations_r100_maf025_multiSNP_phylo.ALL-pars.astral.rooted.tre")

##Run snaq with hmax=0, using astral tree as starting tree, to get starting tree for network inference below..
net0 = snaq!(astraltree_rooted,iqtreeCF_BS10, hmax=0, filename="net0", seed=1234)
#plot(net0, :R); #Sanity Check

##Run snaq with different hmax values.

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net1 = snaq!(net0, iqtreeCF_BS10, hmax=1, filename="net1_bs10", seed=3456)
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net2 = snaq!(net0, iqtreeCF_BS10, hmax=2, filename="net2_bs10", seed=3456)
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net3 = snaq!(net0, iqtreeCF_BS10, hmax=2, filename="net3_bs10", seed=3456)
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net4 = snaq!(net0, iqtreeCF_BS10, hmax=2, filename="net4_bs10", seed=3456)




net1 = readSnaqNetwork("net1_bs10.out") #will read what is generated with snaq! command.
net0 = readSnaqNetwork("net0.out") #read in files after tree generation
net2 = readSnaqNetwork("net2_bs10.out")
net3 = readSnaqNetwork("net3_bs10.out")
net4 = readSnaqNetwork("net4_bs10.out")

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(net0,iqtreeCF_BS10)
df_wide_net0 = fittedQuartetCF(iqtreeCF_BS10)
df_long_net0 = fittedQuartetCF(iqtreeCF_BS10, :long)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(net1,iqtreeCF_BS10) #Rewrites iqtreeCF, so need to run for each network before "fittedQuartetCF" command
df_wide_net1 = fittedQuartetCF(iqtreeCF_BS10)
df_long_net1 = fittedQuartetCF(iqtreeCF_BS10, :long)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(net2,iqtreeCF_BS10)
df_wide_net2 = fittedQuartetCF(iqtreeCF_BS10)
df_long_net2 = fittedQuartetCF(iqtreeCF_BS10, :long)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(net3,iqtreeCF_BS10)
df_wide_net3 = fittedQuartetCF(iqtreeCF_BS10)
df_long_net3 = fittedQuartetCF(iqtreeCF_BS10, :long)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(net4,iqtreeCF_BS10)
df_wide_net4 = fittedQuartetCF(iqtreeCF_BS10)
df_long_net4 = fittedQuartetCF(iqtreeCF_BS10, :long)

CSV.write("fittedCF_bs10_net0.csv",df_long_net0)
CSV.write("fittedCF_bs10_net1.csv",df_long_net1)
CSV.write("fittedCF_bs10_net2.csv",df_long_net2)
CSV.write("fittedCF_bs10_net3.csv",df_long_net3)
CSV.write("fittedCF_bs10_net4.csv",df_long_net4)


###Plot networks, reroot as needed, and determine if networks conflict with outgroup.
netlist1 = readMultiTopology("net1_bs10.networks")
plot(netlist1[1], :R, showGamma=true, showEdgeNumber=true, tipOffset=0.1)
rootonedge!(netlist1[1],17)
plot(netlist1[1], :R, showGamma=true, tipOffset=0.1)

netlist2 = readMultiTopology("net2_bs10.networks")
plot(netlist2[1], :R, showGamma=true, showEdgeNumber=true,tipOffset=0.1)
rootonedge!(netlist2[1],8)
plot(netlist2[1], :R, showGamma=true, tipOffset=0.1, showEdgeNumber=true)

netlist3 = readMultiTopology("net3_bs10.networks")
plot(net3, :R, showGamma=true, showEdgeNumber=true)
rootonedge!(netlist3[1],8)
plot(netlist3[1], :R, showGamma=true, tipOffset=0.1, showEdgeNumber=true)

netlist4 = readMultiTopology("net4_bs10.networks")
plot(net4, :R, showGamma=true, showEdgeNumber=true)
rootonedge!(netlist4[1],8)
plot(netlist4[1], :R, showGamma=true, tipOffset=0.1, showEdgeNumber=true)


##Creating plots with networks that don't conflict with true root position
#Note: not needed with most up to date analyses (i.e., it has not affect because the best networks don't conflict with the root)
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(netlist2[1],iqtreeCF_BS10)
df_wide_net2_secondTree = fittedQuartetCF(iqtreeCF_BS10)
df_long_net2_secondTree = fittedQuartetCF(iqtreeCF_BS10, :long)
CSV.write("fittedCF_bs10_net2_secondTree.csv",df_long_net4)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(netlist3[1],iqtreeCF_BS10)
df_wide_net3_secondTree = fittedQuartetCF(iqtreeCF_BS10)
df_long_net3_secondTree = fittedQuartetCF(iqtreeCF_BS10, :long)
CSV.write("fittedCF_bs10_net3_secondTree.csv",df_long_net4)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(netlist4[1],iqtreeCF_BS10)
df_wide_net4_secondTree = fittedQuartetCF(iqtreeCF_BS10)
df_long_net4_secondTree = fittedQuartetCF(iqtreeCF_BS10, :long)
CSV.write("fittedCF_bs10_net4_secondTree.csv",df_long_net4)


scores = [net0.loglik, net1.loglik, net2.loglik, net3.loglik, net4.loglik] ##For second best networks to fit with outgroups, likelihood added manually
hmax = collect(0:4)
R"plot"(hmax, scores, type="b")
