using PhyloNetworks
#using Pkg
#Pkg.add(url="https://github.com/JuliaPhylo/SNaQ.jl", rev="treechild-galled")
using SNaQ ##This needs to be from the treechild-galled branch 
using DataFrames
using CSV


simulation_dir = "./sims/nets/"
sim_folders = filter(isdir, readdir(simulation_dir, join=true))


for folder in sim_folders
    println(folder)
    net = readnewick(folder*"/net.tre")
    dat_CF = SNaQ.computeexpectedDataCF(net) #compute expected CFs
    df= DataFrame(tablequartetCF(dat_CF)) #convert to a DF 
    CSV.write(folder*"/CFs.csv",df)

    #compute the amount of ILS. Here we compute the deviation of the expected CFs from what we'd expect the CFs to be with no ILS
    #Formally this is average normalized manhatten distance across CFs
    #We 'cheat' on computing the no-ILS case so really it is an approximation
    ILS_deviation = compute_ILS_deviation(net)
    #save result
    open(folder*"/mean_ils_distance.txt", "w") do io
        write(io, string(ILS_deviation))
    end
end

# end

#Not deviation from the No ILS expectation measures the amount of ILS
#This 'cheats' at computing the no ILS expectation by set branches to be long s.t. not coalescing is extremely unlikely
function compute_ILS_deviation(net)

    dat_CF = SNaQ.computeexpectedDataCF(net) #compute expected CFs
    df= DataFrame(tablequartetCF(dat_CF)) #convert to a DF 

    net_no_ILS = deepcopy(net)
    (x->x.length = 100).(net_no_ILS.edge) #set all edges to something that would effectively always coalesce

    dat_CF = SNaQ.computeexpectedDataCF(net_no_ILS) #compute expected CFs under no ILS case
    df_no_ILS= DataFrame(tablequartetCF(dat_CF)) #convert to a DF 
    
    merged_df = innerjoin(df, df_no_ILS, on = [:t1, :t2, :t3, :t4], makeunique = true)

    merged_df.manhattan_dist = abs.(merged_df.CF12_34 .- merged_df.CF12_34_1) .+
                           abs.(merged_df.CF13_24 .- merged_df.CF13_24_1) .+
                           abs.(merged_df.CF14_23 .- merged_df.CF14_23_1)
    merged_df.tv = merged_df.manhattan_dist ./ 2.0

    mean_tv = sum(merged_df.tv)/nrow(merged_df)
end