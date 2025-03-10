point::Int64 = 6601
jobSplit::Int64 = 5
numberOfJobs::Int64 = 1000

for i in 1:numberOfJobs
    run(`qsub -e CH3OH_3D_MEP_AdaptiveGrid_$(point)-$(point + jobSplit).e -o CH3OH_3D_MEP_AdaptiveGrid_$(point)-$(point + jobSplit).o -l h_rt="11:59:00" RunMolproJobs.csh $(point) $(point + jobSplit)`)
    global point = point + jobSplit
    println(point)
end