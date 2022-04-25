using Distributed
addprocs(10)

using Pkg
using PhyloNetworks
using PhyloPlots
using CSV
using RCall
using QuartetNetworkGoodnessFit, DataFrames, CSV
using Gadfly
cd("C:\\Users\\nwhelan\\Documents\\Patera\\RADseq\\SNAQ\\r100-update")

iqtreeTrees_BS10=joinpath("patera_m3M5N5_populations_r100_maf025_multiSNP_nogaps-BS10.trees")
genetrees_BS10=readMultiTopology(iqtreeTrees_BS10)
genetrees[3]  #Sanity check
plot(genetrees[3], :R) #Sainty Check.
q,t = countquartetsintrees(genetrees_BS10) #Generates quartet values (not needed if done with BUCKy)
df = writeTableCF(q,t) #Converts quartets to table.
CSV.write("tableCF_BS10.csv", df) #writes table for use later
iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv") #probably could use df, but wrote to disk and now read back in.
#less("tableCF.csv") #Can be done for sanity check.

##Read in ASTRAL Tree.
astralfile_rooted = ("patera_m3M5N5_populations_r100_maf025_multiSNP_nogaps-rooted.ASTRAL.tre")
astraltree_rooted = readTopology(astralfile)


##Run snaq with hmax=0, using astral tree as starting tree, to get starting tree for network inference below..
net0 = snaq!(astraltree_rooted,iqtreeCF_BS10, hmax=0, filename="net0", seed=1234)
plot(net0, :R); #Sanity Check

##Run snaq with different hmax values.

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net1 = snaq!(net0, iqtreeCF_BS10, hmax=1, filename="net1_bs10", seed=3456)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net2 = snaq!(net0, iqtreeCF_BS10, hmax=2, filename="net2_bs10", seed=3456)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
net3 = snaq!(net0, iqtreeCF_BS10, hmax=2, filename="net3_bs10", seed=3456)

scores = [net0.loglik, net1.loglik, net2.loglik, net3.loglik]
hmax = collect(0:3)
R"plot"(hmax, scores, type="b", ylab="network score", xlab="hmax", col="blue")

###Below lines with read in networks if not in memory.
#net0 = readSnaqNetwork("net0.out") #read in files after tree generation
#net1 = readSnaqNetwork("net1_bs10.out") 
#net2 = readSnaqNetwork("net2_bs10.out")
#net3 = readSnawNetwork(“net3_bs10.out”)


iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv") #Rewrites iqtreeCF, so need to run for each network before "fittedQuartetCF" command
topologyMaxQPseudolik!(net0,iqtreeCF_BS10)
df_wide_net0 = fittedQuartetCF(iqtreeCF_BS10)
df_long_net0 = fittedQuartetCF(iqtreeCF_BS10, :long)

iqtreeCF_BS10 = readTableCF("tableCF_BS10.csv")
topologyMaxQPseudolik!(net1,iqtreeCF_BS10)
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


##See R code for plotting. Plotting withing Julia kept failing on my system.
CSV.write("fittedCF_bs10_net0.csv",df_long_net0)
CSV.write("fittedCF_bs10_net1.csv",df_long_net1)
CSV.write("fittedCF_bs10_net2.csv",df_long_net2)
CSV.write("fittedCF_bs10_net3.csv",df_long_net3)

###Plotting commands for phylogenetic networks
plot(net1, :R, showGamma=true, showEdgeNumber=true);
rootonedge!(net1,6);
plot(net2, :R, showGamma=true, showEdgeNumber=true);