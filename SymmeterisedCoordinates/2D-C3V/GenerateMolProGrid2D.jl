using Printf

stretchesSpacing::Vector{Float64} = [0.0000, -0.025, 0.025, 0.050, -0.0500, -0.1000, 0.1000, 0.12500, -0.125000, 0.1500, -0.1500, 0.200, 
   -0.200, 0.300, -0.300, 0.400, 0.500, 0.600, 0.700]
angleSpacing::Vector{Float64} = [0.0000, 1.000, -1.000, 2.500, -2.5000, 5.0000, -5.0000, 7.500, -7.500, 10.0000, -10.0000, 20.0000, 
    -20.0000, 30.0000, -30.0000, 50.0000, -50.0000, 60.0000]
dihedralSpacing::Vector{Float64} = [0.0000, 1.000, -1.000, 2.5000, -2.5000, 5.0000, -5.0000, 7.500, -7.500, 10.0000, -10.0000, 40.0000, 
    -40.0000, 60.0000, -60.0000, 80.0000, -80.0000, 100.0000]
torsionSpacing::Vector{Float64} = [0.0000, 5.0000, 10.0000, 15.0000, 20.00000, 25.000, 30.0000, 35.0000, 40.0000, 45.0000, 50.0000, 55.0000, 60.0000]

stretchesGrid::Int64 = size(stretchesSpacing)[1]
angleGrid::Int64 = size(angleSpacing)[1]
dihedralGrid::Int64 = size(angleSpacing)[1]
torsionGrid::Int64 = size(torsionSpacing)[1]
convertToRadians::Float64 = 2*pi/360

function EqCH(tau::Float64)::Float64
    tau = tau*convertToRadians
    a0::Float64 = 1.08902324e+00
    a1::Float64 = 2.71787366e-03
    a2::Float64 = -2.16000169e-03
    a3::Float64 = -2.55337930e-04
    rCH::Float64 = a0 + a1*cos(tau) + a2*cos(2*tau) + a3*cos(3*tau)
    return 1.08924348 #rCH
end

function EqaHCO(tau::Float64)::Float64
    tau = tau*convertToRadians
    a0::Float64 = 1.10243240e+02 
    a1::Float64 = 2.29302542e+00
    a2::Float64 = -7.16664340e-01
    a3::Float64 = 7.93169192e-02
    a4::Float64 = 3.31054618e-02
    rHCO::Float64 = a0 + a1*cos(tau) + a2*cos(2*tau) + a3*cos(3*tau) + a4*cos(4*tau)
    return 110.16192385 # rHCO
end

function ConvertFromSymmeterisedStretches(stretchesDisplacement::Vector{Float64}, equilibriumBondLengths::Vector{Float64})::Vector{Float64}
    stretches::Vector{Float64} = equilibriumBondLengths
    stretches[1] += sqrt(2/3)*stretchesDisplacement[1] + stretchesDisplacement[2]/sqrt(6) + stretchesDisplacement[3]/sqrt(2)
    stretches[2] += sqrt(2/3)*stretchesDisplacement[1] - sqrt(2/3)*stretchesDisplacement[2]
    stretches[3] += sqrt(2/3)*stretchesDisplacement[1] + stretchesDisplacement[2]/sqrt(6) - stretchesDisplacement[3]/sqrt(2)
    return stretches
end

function ConvertFromSymmeterisedAngles(anglesDisplacement::Vector{Float64}, equilibriumAngles::Vector{Float64})::Vector{Float64}
    angles::Vector{Float64} = equilibriumAngles
    angles[1] += sqrt(2/3)*anglesDisplacement[1] - anglesDisplacement[2]/sqrt(6) - anglesDisplacement[3]/sqrt(2)
    angles[2] += sqrt(2/3)*anglesDisplacement[1] + sqrt(2/3)*anglesDisplacement[2]
    angles[3] += sqrt(2/3)*anglesDisplacement[1] - anglesDisplacement[2]/sqrt(6) + anglesDisplacement[3]/sqrt(2)
    return angles
end

d1::Float64 =  60.0  # 61.43364279
d2::Float64 =  180.0 #180.00000000
d3::Float64 =  300.0 #298.56635721
SaEq::Float64 = 0 # -1.7558466544600673
SbEq::Float64 = 0 #  3.04121561582463
tau::Float64 = 60.00

ahh1::Float64 = tau+1.0/3.0*sqrt(2.0)*SbEq
ahh2::Float64 = 120.0+tau-1.0/6.0*sqrt(2.0)*SbEq-1.0/6.0*sqrt(6.0)*SaEq
ahh3::Float64 = 240.0+tau-1.0/6.0*sqrt(2.0)*SbEq+1.0/6.0*sqrt(6.0)*SaEq

rCOeq::Float64 =                1.42077677
rOHeq::Float64=                 0.96013932
aCOHeq::Float64=              108.12930637
aHH1eq::Float64=              60.0 # 61.43364279
aHH2eq::Float64=              180.0 #180.00000000
aHH3eq::Float64=              300.0 #298.56635721
tauEq::Float64 = 60.00000
rCH1eq::Float64=             EqCH(tauEq)         #1.09108970
rCH2eq::Float64=             EqCH(tauEq + 120)   #1.08555104
rCH3eq::Float64=             EqCH(tauEq + 240)    #1.09108970
aOCH1eq::Float64=            EqaHCO(tauEq)        #111.95221297
aOCH2eq::Float64=            EqaHCO(tauEq + 120)  #106.58134561
aOCH3eq::Float64=            EqaHCO(tauEq + 240)  #111.95221297
symmeterisedStretchEq::Float64 = 0.00
symmeterisedStretchAEq::Float64 = 0.00
symmeterisedStretchBEq::Float64 = 0.00
symmeterisedAnglesEq::Float64 = 0.00
symmeterisedAnglesAEq::Float64 = 0.00
symmeterisedAnglesBEq::Float64 = 0.00
equilibriumGrid::Vector{Float64} = [rCOeq, rOHeq, rCH1eq, rCH2eq, rCH3eq, aCOHeq, aOCH1eq, aOCH2eq, aOCH3eq, aHH1eq, aHH2eq, aHH3eq]
equilibriumSymmeterisedGrid::Vector{Float64} = [rCOeq, rOHeq, r1SymmeterisedEq, r2SymmeterisedEq, r3SymmeterisedEq, aCOHeq, a1SymmeterisedEq, a2SymmeterisedEq, a3SymmeterisedEq, SaEq, SbEq, tauEq]

function PrintGeometry(point::Int64, grid::Vector{Float64})
    @printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", point, grid[1], grid[2], grid[3], grid[4], grid[5], grid[6], grid[7], grid[8], grid[9], grid[10], grid[11], grid[12])
end

function SubmitJob(point::Int64, grid::Vector{Float64})
    submission::Cmd = `qsub -e CH3OH_2D_MEP_AdaptiveGrid_$(point).e -o CH3OH_2D_MEP_AdaptiveGrid_$(point).o -l h_rt="11:59:00" GenerateMolproScript2D.csh $(point) $(grid[1]) $(grid[2]) $(grid[3]) $(grid[4]) $(grid[5]) $(grid[6]) $(grid[7]) $(grid[8]) $(grid[9]) $(grid[10]) $(grid[11]) $(grid[12])`
    run(submission)
end

point::Int64 = 206
grids::Vector{Vector{Float64}} = []

for i in 1:12
    for j in i+1:12
        for k in 1:size(spacingOfGrid[i])[1]
            for l in 1:size(spacingOfGrid[j])[1]
                displacementVector::Vector{Float64} = zeros(12)
                displacementVector[i] = spacingOfGrid[i][k]
                displacementVector[j] = spacingOfGrid[j][l]
                gridSymmeterised::Vector{Float64} = equilibriumSymmeterisedGrid + displacementVector
                grid::Vector{Float64} = gridSymmeterised .+ 0
                equilibriumBondLengths::Vector{Float64} = [EqCH(gridSymmeterised[12]), EqCH(gridSymmeterised[12]+120), EqCH(gridSymmeterised[12]+240)]
                equilibriumBondAngles::Vector{Float64} = [EqaHCO(gridSymmeterised[12]), EqaHCO(gridSymmeterised[12]+120), EqaHCO(gridSymmeterised[12]+240)]
                # println("$(i) $(j) $(displacementVector[i]) $(displacementVector[j])")
                # println(gridSymmeterised[10:11])
                grid[3:5] = ConvertFromSymmeterisedStretches(grid[3:5], equilibriumBondLengths)
                grid[7:9] = ConvertFromSymmeterisedAngles(grid[7:9], equilibriumBondAngles)
                grid[10] = gridSymmeterised[12]+1.0/3.0*sqrt(2.0)*gridSymmeterised[11]
                grid[11] = 120.0+gridSymmeterised[12]-1.0/6.0*sqrt(2.0)*gridSymmeterised[11]-1.0/6.0*sqrt(6.0)*gridSymmeterised[10]
                grid[12] = 240.0+gridSymmeterised[12]-1.0/6.0*sqrt(2.0)*gridSymmeterised[11]+1.0/6.0*sqrt(6.0)*gridSymmeterised[10]
                push!(grids, grid)
                PrintGeometry(point, grid)
                global point = point + 1
            end
        end
    end
end

# for i in 1:5
#     for k in i+1:5
#         for j in 2:stretchesGrid
#             displacementVector::Vector{Float64} = zeros(12)
#             displacementVector[i] = stretchesSpacing[j]
#             for l in 2:stretchesGrid
#                 displacementVector[k] = stretchesSpacing[l]
#                 grid::Vector{Float64} = equilibriumGrid + displacementVector
#                 push!(grids, grid)
#                 PrintGeometry(point, grid)
#                 global point = point + 1
#             end
#         end
#     end
#     for k in 6:9
#         for j in 2:stretchesGrid
#             displacementVector::Vector{Float64} = zeros(12)
#             displacementVector[i] = stretchesSpacing[j]
#             for l in 2:angleGrid
#                 displacementVector[k] = angleSpacing[l]
#                 grid::Vector{Float64} = equilibriumGrid + displacementVector
#                 push!(grids, grid)
#                 PrintGeometry(point, grid)
#                 global point = point + 1
#             end
#         end
#     end
#     for k in 10:11
#         for j in 2:stretchesGrid
#             displacementVector::Vector{Float64} = zeros(12)
#             displacementVector[i] = stretchesSpacing[j]
#             for l in 2:dihedralGrid
#                 symmeterisedDihedrals::Vector{Float64} = [-1.7558466544600673, 3.04121561582463]
#                 symmeterisedDihedrals[k - 9] = symmeterisedDihedrals[k - 9] + dihedralSpacing[l]
#                 grid::Vector{Float64} = equilibriumGrid + displacementVector
#                 grid[10] = tauEq+1.0/3.0*sqrt(2.0)*symmeterisedDihedrals[2]
#                 grid[11] = 120.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]-1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
#                 grid[12] = 240.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]+1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
#                 PrintGeometry(point, grid)
#                 push!(grids, grid)
#                 global point = point + 1
#             end
#         end
#     end
#     for j in 2:stretchesGrid
#         displacementVector::Vector{Float64} = zeros(12)
#         displacementVector[i] = stretchesSpacing[j]
#         for l in 2:torsionGrid
#             displacementVector[10] = torsionSpacing[l]
#             displacementVector[11] = torsionSpacing[l]
#             displacementVector[12] = torsionSpacing[l]
#             grid::Vector{Float64} = equilibriumGrid + displacementVector
#             grid[3] = EqCH(tauEq + torsionSpacing[l])
#             grid[4] = EqCH(tauEq + torsionSpacing[l] + 120)
#             grid[5] = EqCH(tauEq + torsionSpacing[l] + 240)
#             grid[7] = EqaHCO(tauEq + torsionSpacing[l])
#             grid[8] = EqaHCO(tauEq + torsionSpacing[l] + 120)
#             grid[9] = EqaHCO(tauEq + torsionSpacing[l] + 240)
#             PrintGeometry(point, grid)
#             push!(grids, grid)
#             global point = point + 1
#         end
#     end
# end

# for i in 6:9
#     for k in i+1:9
#         for j in 2:angleGrid
#             displacementVector::Vector{Float64} = zeros(12)
#             displacementVector[i] = angleSpacing[j]
#             for l in 2:angleGrid
#                 displacementVector[k] = angleSpacing[l]
#                 grid::Vector{Float64} = equilibriumGrid + displacementVector
#                 push!(grids, grid)
#                 PrintGeometry(point, grid)
#                 global point = point + 1
#             end
#         end
#     end
#     for k in 10:11
#         for j in 2:angleGrid
#             displacementVector::Vector{Float64} = zeros(12)
#             displacementVector[i] = angleSpacing[j]
#             for l in 2:dihedralGrid
#                 symmeterisedDihedrals::Vector{Float64} = [-1.7558466544600673, 3.04121561582463]
#                 symmeterisedDihedrals[k - 9] = symmeterisedDihedrals[k - 9] + dihedralSpacing[l]
#                 grid::Vector{Float64} = equilibriumGrid + displacementVector
#                 grid[10] = tauEq+1.0/3.0*sqrt(2.0)*symmeterisedDihedrals[2]
#                 grid[11] = 120.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]-1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
#                 grid[12] = 240.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]+1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
#                 PrintGeometry(point, grid)
#                 push!(grids, grid)
#                 global point = point + 1
#             end
#         end
#     end
#     for j in 2:angleGrid
#         displacementVector::Vector{Float64} = zeros(12)
#         displacementVector[i] = angleSpacing[j]
#         for l in 2:torsionGrid
#             displacementVector[10] = torsionSpacing[l]
#             displacementVector[11] = torsionSpacing[l]
#             displacementVector[12] = torsionSpacing[l]
#             grid::Vector{Float64} = equilibriumGrid + displacementVector
#             grid[3] = EqCH(tauEq + torsionSpacing[l])
#             grid[4] = EqCH(tauEq + torsionSpacing[l] + 120)
#             grid[5] = EqCH(tauEq + torsionSpacing[l] + 240)
#             grid[7] = EqaHCO(tauEq + torsionSpacing[l])
#             grid[8] = EqaHCO(tauEq + torsionSpacing[l] + 120)
#             grid[9] = EqaHCO(tauEq + torsionSpacing[l] + 240)
#             PrintGeometry(point, grid)
#             push!(grids, grid)
#             global point = point + 1
#         end
#     end
# end

# for i in 10:11
#     for k in i+1:11
#         for j in 2:dihedralGrid
#             displacementVector::Vector{Float64} = zeros(12)
#             symmeterisedDihedrals::Vector{Float64} = [-1.7558466544600673, 3.04121561582463]
#             symmeterisedDihedrals[i - 9] = symmeterisedDihedrals[i - 9] + dihedralSpacing[j]
#             for l in 2:dihedralGrid
#                 symmeterisedDihedrals[k - 9] = symmeterisedDihedrals[k - 9] + dihedralSpacing[l]
#                 grid::Vector{Float64} = equilibriumGrid + displacementVector
#                 grid[10] = tauEq+1.0/3.0*sqrt(2.0)*symmeterisedDihedrals[2]
#                 grid[11] = 120.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]-1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
#                 grid[12] = 240.0+tauEq-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]+1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
#                 PrintGeometry(point, grid)
#                 push!(grids, grid)
#                 global point = point + 1
#             end
#         end
#     end
#     for j in 2:dihedralGrid
#         displacementVector::Vector{Float64} = zeros(12)
#         symmeterisedDihedrals::Vector{Float64} = [-1.7558466544600673, 3.04121561582463]
#         symmeterisedDihedrals[i - 9] = symmeterisedDihedrals[i - 9] + dihedralSpacing[j]
#         for l in 2:torsionGrid
#             grid::Vector{Float64} = equilibriumGrid + displacementVector
#             grid[3] = EqCH(tauEq + torsionSpacing[l])
#             grid[4] = EqCH(tauEq + torsionSpacing[l] + 120)
#             grid[5] = EqCH(tauEq + torsionSpacing[l] + 240)
#             grid[7] = EqaHCO(tauEq + torsionSpacing[l])
#             grid[8] = EqaHCO(tauEq + torsionSpacing[l] + 120)
#             grid[9] = EqaHCO(tauEq + torsionSpacing[l] + 240)
#             grid[10] = tauEq+torsionSpacing[l]+1.0/3.0*sqrt(2.0)*symmeterisedDihedrals[2]
#             grid[11] = 120.0+tauEq+torsionSpacing[l]-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]-1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
#             grid[12] = 240.0+tauEq+torsionSpacing[l]-1.0/6.0*sqrt(2.0)*symmeterisedDihedrals[2]+1.0/6.0*sqrt(6.0)*symmeterisedDihedrals[1]
#             PrintGeometry(point, grid)
#             push!(grids, grid)
#             global point = point + 1
#         end
#     end
# end