point::Int64 = 18081
jobSplit::Int64 = 5
numberOfJobs::Int64 = 1000

for i in 1:numberOfJobs
    run(`qsub -e CH3OH_CBS_$(point)-$(point + jobSplit).e -o CH3OH_CBS_$(point)-$(point + jobSplit).o -l h_rt="5:59:00" -l mem=70G -l tmpfs=100G RunMolproJobs.csh $(point) $(point + jobSplit)`)
    global point = point + jobSplit
    println(point)
end