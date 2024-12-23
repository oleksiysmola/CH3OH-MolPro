point::Int64 = 18415
jobSplit::Int64 = 6
numberOfJobs::Int64 = 1000

for i in 1:numberOfJobs
    run(`qsub -e CH3OH_2D_MEP_AdaptiveGrid_$(point)-$(point + jobSplit).e -o CH3OH_2D_MEP_AdaptiveGrid_$(point)-$(point + jobSplit).o -l h_rt="11:59:00" -l mem=30G -l tmpfs=100G RunMolproJobs.csh $(point) $(point + jobSplit)`)
    global point = point + jobSplit
    println(point)
end