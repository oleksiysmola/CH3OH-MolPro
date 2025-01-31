#!/bin/csh
#
# Generation of the input file for MOLPRO
#

set pwd = `pwd`

set point = $1        
set directory = /scratch/scratch/zcaposm/Methanol/MOLPRO/ExtrapolationCBS/CH3OH_CBS
set fname = CH3OH_CBS_${point}

cat<<endb> ${fname}.inp
***, Methanol ground state at CCSD(T)-F12 and aug-cc-pV(T-Q)Z-F12 to CBS
memory,1000,m;
gthresh,energy=1.d-10,zero=1.d-14,thrint=1.d-14,oneint=1.d-14,twoint=1.d-14,prefac=1.d-20

geometry={angstrom
c 
o , 1, rco 
h , 2, roh, 1, acoh
h1, 1, rch1, 2, aoch1  3  ahh1
h2, 1, rch2, 2, aoch2, 3, ahh2
h3, 1, rch3, 2, aoch3, 3, ahh3
}

rco= $2
roh= $3
rch1= $4
rch2= $5
rch3= $6
acoh= $7
aoch1= $8
aoch2= $9
aoch3= $10
ahh1=$11
ahh2=$12
ahh3=$13

basis=aug-cc-pVTZ-f12
hf
frozen
ccsd(t)-f12,ri_basis=cc-pV5Z/jkfit,df_basis=aug-cc-pwCV5Z/mp2fit

E_ccsd3 = ENERGC(2)
E_t3 = ENERGY(2) - ENERGC(2)

basis=aug-cc-pVQZ-f12
hf
frozen
ccsd(t)-f12,ri_basis=cc-pV5Z/jkfit,df_basis=aug-cc-pwCV5Z/mp2fit

E_ccsd4 = ENERGC(2)
E_t4 = ENERGY(2) - ENERGC(2)

E_ccsd = (E_ccsd4 - E_ccsd3)*1.363388 + E_ccsd3
E_t = (E_t4 - E_t3)*1.769474 + E_t3
E_total = E_ccsd + E_t

! Output the energy
xxx = "mmm"
point = ${point}
text ### CH3OH
table,xxx,rco,roh,rch1,rch2,rch3,acoh,aoch1,aoch2,aoch3,ahh1,ahh2,ahh3,E_total,point
DIGITS, 0, 8, 8, 8, 8, 8, 8, 8, 8, 8,  8,  8,  8, 12, 4

--- End of Script ---
endb

# module load molpro/2020.1/openmp
molpro ${fname}.inp
rm ${fname}.inp
cp ${fname}.out ${directory}
