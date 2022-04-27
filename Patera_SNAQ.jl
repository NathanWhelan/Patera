using Distributed
addprocs(10)

using Pkg
using PhyloNetworks
using PhyloPlots
using RCall
using QuartetNetworkGoodnessFit
using DataFrames
using CSV
using Gadfly
cd("C:\\Users\\nwhelan\\Documents\\Patera\\RADseq\\SNAQ\\r100-update")

iqtreeTrees=joinpath("patera_m3M5N5_populations_r100_maf025_multiSNP_nogaps.trees")
genetrees=readMultiTopology(iqtreeTrees)
#genetrees[3]  #Sanity check
#plot(genetrees[3], :R) #Sainty Check.
q,t = countquartetsintrees(genetrees) #Generates quartet values (not needed if done with BUCKy)
df = writeTableCF(q,t) #Converts quartets to table.
CSV.write("tableCF.csv", df) #writes table for use later
iqtreeCF = readTableCF("tableCF.csv") #probably could use df, but wrote to disk and now read back in.
#less("tableCF.csv") #Can be done for sanity check.

iqtreeTrees_BS10=joinpath("patera_m3M5N5_populations_r100_maf025_multiSNP_nogaps-BS10.trees")
genetrees_BS10=readMultiTopology(iqtreeTrees_BS10)
#genetrees[3]  #Sanity check
#plot(genetrees[3], :R) #Sainty Check.
q,t = countquartetsintrees(genetrees_BS10) #Generates quartet values (not needed if done with BUCKy)
df = writeTableCF(q,t) #Converts quartets to table.
CSV.write("tableCF_BS10.csv", df) #writes table for use later
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv") #probably could use df, but wrote to disk and now read back in.
#less("tableCF.csv") #Can be done for sanity check.

##Read in ASTRAL tree for starting tree.
astralfile_rooted = ("patera_m3M5N5_populations_r100_maf025_multiSNP_nogaps-rooted.ASTRAL.tre")
astraltree_rooted = readTopology(astralfile_rooted)


##Run snaq with hmax=0, using astral tree as starting tree, to get starting tree for network inference below..
net0 = snaq!(astraltree_rooted,iqtreeCF_BS10, hmax=0, filename="net0", seed=1234)
plot(net0, :R); #Sanity Check

##Run snaq with different hmax values.

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net1_test = snaq!(net0, iqtreeCF_BS10, hmax=1, filename="net1_bs10", seed=3456)
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net2_test = snaq!(net0, iqtreeCF_BS10, hmax=2, filename="net2_bs10", seed=3456)
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net3_test = snaq!(net0, iqtreeCF_BS10, hmax=2, filename="net3_bs10", seed=3456)
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net4_test = snaq!(net0, iqtreeCF_BS10, hmax=2, filename="net4_bs1", seed=3456)



net0 = readSnaqNetwork("net0.out") #read in files after tree generation
net1 = readSnaqNetwork("net1_bs10.out") #will read what is generated with snaq! command.
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
##See R code for plotting expected vs observed CF. Plotting withing Julia kept failing on my system.

###Plot networks, reroot as needed, and determine if networks conflict with outgroup.
netlist1 = readMultiTopology("net1.networks")
plot(netlist1[1], :R, showGamma=true, showEdgeNumber=true, tipOffset=0.1)
rootonedge!(netlist1[1],7)
plot(netlist1[1], :R, showGamma=true, tipOffset=0.1)

netlist2 = readMultiTopology("net2.networks")
plot(netlist2[2], :R, showGamma=true, tipOffset=0.1)

plot(net3, :R, showGamma=true, showEdgeNumber=true)
netlist3 = readMultiTopology("net3.networks")
rootonedge!(netlist3[2],8)
plot(netlist3[2], :R, showGamma=true, tipOffset=0.1, showEdgeNumber=true)

plot(net4, :R, showGamma=true, showEdgeNumber=true)
netlist4 = readMultiTopology("net3.networks")
rootonedge!(netlist4[2],8)
plot(netlist4[2], :R, showGamma=true, tipOffset=0.1, showEdgeNumber=true)


##Creating plots with networks that don't conflict with true root position
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(netlist2[2],iqtreeCF_BS10)
df_wide_net2_secondTree = fittedQuartetCF(iqtreeCF_BS10)
df_long_net2_secondTree = fittedQuartetCF(iqtreeCF_BS10, :long)
CSV.write("fittedCF_bs10_net2_secondTree.csv",df_long_net4)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(netlist3[2],iqtreeCF_BS10)
df_wide_net3_secondTree = fittedQuartetCF(iqtreeCF_BS10)
df_long_net3_secondTree = fittedQuartetCF(iqtreeCF_BS10, :long)
CSV.write("fittedCF_bs10_net3_secondTree.csv",df_long_net4)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(netlist4[2],iqtreeCF_BS10)
df_wide_net4_secondTree = fittedQuartetCF(iqtreeCF_BS10)
df_long_net4_secondTree = fittedQuartetCF(iqtreeCF_BS10, :long)
CSV.write("fittedCF_bs10_net4_secondTree.csv",df_long_net4)


scores = [net0.loglik, net1.loglik, 25.65196551371749, 25.65196551371749, 25.65196551371749] ##For second best networks to fit with outgroups, likelihood added manually
hmax = collect(0:4)
R"plot"(hmax, scores, type="b")



