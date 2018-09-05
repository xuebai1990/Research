#! /bin/bash -l
#PBS -l nodes=1:ppn=2,walltime=120:00:00
#PBS -N C4S_1
#PBS -q E5v4

# For metropolis
# Use melting file as input file
# Box 1 should be vapor box
#*************************************************************************************
# Change these values:
  RUN=1
  SEED=1001
  TEMP=280
  NMOL=240      # number of molecules
  STEPC=5000    # number of steps for cooling
  STEPE=50000   # number of steps for equilibrium
  STEPP=100000  # number of steps for production
  MPI=2         # number of mpi run
#*************************************************************************************

DATA=$PBS_O_WORKDIR
cd $DATA
FILELOG=log

PROG=$HOME/exe-metropolis2-precision/src/topmon

# Change seed
OLD1="seed="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "seed"){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
NEW2="$SEED"
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4

# Melting
echo "Melting: " >> $FILELOG
/usr/bin/time -ao $FILELOG $PROG || exit -1
cp -f config1a.dat fort.77
mv -f run1a.dat run1a.melting
mv -f config1a.dat config.melting
mv -f movie1a.dat movie.melting
mv -f fort.12 fort12.melting
for i in 1 2;do
   mv -f "box${i}config1a.xyz" "box${i}config.melting"
   mv -f "box${i}movie1a.xyz" "box${i}movie.melting"
done

# Cooling
echo "Cooling: " >> $FILELOG

# Change linit to false
OLD1="linit="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "linit"){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
NEW2=F
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4

# Change nstep
OLD1="nstep="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "nstep"){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
NEW2="$STEPC"
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4

# Change temperature
STRING1=$(awk '{if($0 ~ "boxlx"){getline;print $0;exit}}' fort.4)
TEMPOLD=$(awk '{if($0 ~ "boxlx"){getline;print $12;exit}}' fort.4)
STRING2=$(echo ${STRING1} | sed "s/${TEMPOLD}/${TEMP}/")
sed -i "s/${STRING1}/${STRING2}/" fort.4
STRING1=$(awk '{if($0 ~ "boxlx"){for(i=1;i<8;i++)getline;print $0;exit}}' fort.4)
TEMPOLD=$(awk '{if($0 ~ "boxlx"){for(i=1;i<8;i++)getline;print $12;exit}}' fort.4)
STRING2=$(echo ${STRING1} | sed "s/${TEMPOLD}/${TEMP}/")
sed -i "s/${STRING1}/${STRING2}/" fort.4


/usr/bin/time -ao $FILELOG mpirun -np "${MPI}" $PROG || exit -1
cp -f config1a.dat fort.77
mv -f run1a.dat run1a.cooling
mv -f config1a.dat config.cooling
mv -f movie1a.dat movie.cooling
mv -f fort.12 fort12.cooling
for i in 1 2;do
   mv -f "box${i}config1a.xyz" "box${i}config.cooling"
   mv -f "box${i}movie1a.xyz" "box${i}movie.cooling"
done

# Change rcut
RCUT=14.0   #change your rcut here
STRING1=$(awk '{if($0 ~ "boxlx"){for(i=1;i<8;i++)getline;print $0;exit}}' fort.4)
RCUTOLD=$(awk '{if($0 ~ "boxlx"){for(i=1;i<8;i++)getline;print $4;exit}}' fort.4)
STRING2=$(echo ${STRING1} | sed "s/${RCUTOLD}/${RCUT}/")
sed -i "s/${STRING1}/${STRING2}/" fort.4

# Set up volume move
PMV=$(echo "2.5/${NMOL}" | bc -l)
OLD1="pmvol="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "pmvol="){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
NEW2="$PMV"
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4

# Change rmin
OLD1="rmin="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "rmin"){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
NEW2=1.20
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4

echo "Equilibrium${i}: " >> $FILELOG

# Change nstep
OLD1="nstep="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "nstep"){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
NEW2="$STEPE"
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4

OLD1="pmswap="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "pmswap"){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
SP=0.4
NEW2=$(echo "${SP}+${PMV}" | bc -l)
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4

/usr/bin/time -ao $FILELOG mpirun -np "${MPI}" $PROG || exit -1
cp -f config1a.dat fort.77
mv -f run1a.dat run1a.equil
mv -f config1a.dat config.equil
mv -f movie1a.dat movie.equil
mv -f fort.12 fort12.equil
for j in 1 2;do
   mv -f "box${j}config1a.xyz" "box${j}config.equil"
   mv -f "box${j}movie1a.xyz" "box${j}movie.equil"
done

# Production

echo "Production: " >> $FILELOG

# Change nstep
OLD1="nstep="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "nstep"){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
NEW2="$STEPP"
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4

# Change iratio iratv
OLD1="iratio="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "iratio"){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
NEW2="$STEPP"
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4
OLD1="iratv="
OLD2=$(awk 'BEGIN {FS="="}{if($0 ~ "iratv"){print $2}}' fort.4)
OLD="${OLD1}${OLD2}"
NEW2="$STEPP"
NEW="${OLD1}${NEW2}"
sed -i "s/${OLD}/${NEW}/" fort.4

/usr/bin/time -ao $FILELOG mpirun -np "${MPI}" $PROG || exit -1
cp -f config1a.dat fort.77
mv -f run1a.dat run1a.prod
mv -f config1a.dat config.prod
mv -f movie1a.dat movie.prod
mv -f fort.12 fort12.prod
for i in 1 2 3;do
   mv -f "box${i}config1a.xyz" "box${i}config.prod"
   mv -f "box${i}movie1a.xyz" "box${i}movie.prod"
done

echo "Finished!" >> $FILELOG


exit 0



