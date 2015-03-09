#! /bin/bash

SIM_HEADER="pre_"

PARTICLE_TYPE='Proton'
EXP_FIT_SOFTWARE=$HOME/bin/exponential_fit
SPEC_DECODER=$HOME/bin/leggi_diag
#il tool seguente e` scan-columns in tools-Propaga
#serve per splittare i risultati collezionati dallo script su diversi files,
# in funzione di un parametro, per produrre piu` plots
SCANNER=$HOME/bin/scan-columns
DO_SCAN=true
#per il seguente, copiarsi dal prepare_scan_input la riga che genera tutti i valori (in questo caso di bulk lengths scannerizzate)
columns_values=$(awk 'BEGIN{for(i=2.0;i<=4.0;i+=1.0)print i}')
column=4
# 1:{PREPLASMA_LENGTH} 2:{DENSITY} 3:{RAMP_LENGTH} 4:{BULK_LENGTH} 5:{CONT_LENGTH} 
# 6:{PROTON_MAX_ENERGY} 7:{PROTON_TOT_ENERGY} 8:{PROTON_AVE_ENERGY} 9:{PROTON_TOT_NUMBER}" 
# 10:{SEL_PROTON_AVE_ENERGY} 11:{SEL_PROTON_TOT_NUMBER}

FINAL_LINE_NUMBER=146
DIAG_STEP_TO_BE_READ=10
SPEC_TIME_TO_BE_READ=100

OUTPUT_FILE="energy_scan.txt"




SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))

rm -f ${OUTPUT_FILE}
touch ${OUTPUT_FILE}


printf "# PreplasmaLength \t Density \t RampLength \t BulkLength \t ${PARTICLE_TYPE}MaxEnergy \t ${PARTICLE_TYPE}TotEnergy \t ${PARTICLE_TYPE}AveEnergy \t ${PARTICLE_TYPE}TotNumber \t Sel${PARTICLE_TYPE}AveEnergy \t Sel${PARTICLE_TYPE}TotNumber \n" >> ${OUTPUT_FILE} 

for sim in "${SIMULATION_FOLDERS[@]}"
do
#pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam}
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 CONT_LENGTH=($(echo $sim | awk -F'_' '{print $10}'))

cd $sim
if [ -f "diag${DIAG_STEP_TO_BE_READ}.dat" ];
then
 PROTON_MAX_ENERGY=($(head -${FINAL_LINE_NUMBER} diag${DIAG_STEP_TO_BE_READ}.dat |tail -1 |awk '{print $4}'))
 PROTON_TOT_ENERGY=($(head -${FINAL_LINE_NUMBER} diag${DIAG_STEP_TO_BE_READ}.dat |tail -1 |awk '{print $3}'))
else
 PROTON_MAX_ENERGY="0.0"
 PROTON_TOT_ENERGY="0.0"
fi


if [ -f "spec${DIAG_STEP_TO_BE_READ}.dat" ];
then
 ${SPEC_DECODER} spec${DIAG_STEP_TO_BE_READ}.dat
 aveData=($( ${EXP_FIT_SOFTWARE} -scan spec${DIAG_STEP_TO_BE_READ}.dat_${PARTICLE_TYPE}_${SPEC_TIME_TO_BE_READ}.txt  ))
else
 aveData[0]="0.0"
 aveData[1]="0.0"
 aveData[2]="0.0"
 aveData[3]="0.0"
fi

PROTON_AVE_ENERGY=${aveData[0]}
PROTON_TOT_NUMBER=${aveData[1]}
SEL_PROTON_AVE_ENERGY=${aveData[2]}
SEL_PROTON_TOT_NUMBER=${aveData[3]}

cd ..

#read -p "Press [Enter] key to continue..."
printf '\t      %.1f \t%.1f\t     %.1f   \t%.1f\t %.1f\t%.2f\t%.3e \t%s\t%s\t%s\t%s' "${PREPLASMA_LENGTH}" "${DENSITY}" "${RAMP_LENGTH}" "${BULK_LENGTH}" "${CONT_LENGTH}" "${PROTON_MAX_ENERGY}" "${PROTON_TOT_ENERGY}" "${PROTON_AVE_ENERGY}" "${PROTON_TOT_NUMBER}" "${SEL_PROTON_AVE_ENERGY}" "${SEL_PROTON_TOT_NUMBER}" >> ${OUTPUT_FILE}

printf '\n'  >> ${OUTPUT_FILE}


done

if [ ${DO_SCAN} ] 
then
 for value in ${columns_values}
 do
  $SCANNER -in energy_scan.txt -out energy_scan_$value.txt -select $column $value
 done
else 
 echo "output on ${OUTPUT_FILE}"
fi