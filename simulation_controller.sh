#! /bin/bash
## MEMBER OF CYCLE 1
nA=3
## MEMEBER OF CYCLE 2
nB=3
## kSetGeneration: 0: does not generate k_set or k_SET; 1: generates k_set; 2: generates k_SET
kSetGeneration=1
## DominantCycle: <0: no catalitic aid,  1: A is the dominant cycle; 2: B is the dominant cycle; 3: A and B cycles coexist; 8: system died out;, 9: non identifiable simulation
DominantCycle=-1
## RANSOM SEED INITIALIZE
InitRandNum=0
## 0: no reverse reaction, 1: one reverse reaction
reverse=0
## ratio: <0: no revers raction, >0: revers rate = toward rate/ratio: (10 9 8 7 6 5 4 3 2 1 0.5)
ratio=(10 9 8 7 6 5 4 3 2 1 0.5)
## HetCat: 0: no cross-catalysis; 1: CC1: sub1: A1 or B1, sub2: X; 2: CC2: sub1: X, sub2: A1 or B1 (CC: monomolecular heterocatalysis selected from the  results of CC1)
HetCat=0
## PozRev: <0: no reverz, >0: pozition of revers reaction (PozRev-1 in the *.c file)
PozRev=-1
## Tmodel: 0: unlimited resources, 1: chemostate
TModel=1
## JoinedCycles 0: no joined cycles, 1: cycles joined at JoinedPointA=JoinedPointB
JoinedCycles=0
## JoinedPointA: 0: no joined cycles, >0: cycle A is joined to cycle B at this point
JoinedPointA=0
## JoinedPointB: 0: no joined cycles, >0: cycle B is joined to cycle A at this point
JoinedPointB=0

## DELETE OUTPUT FILES (previous simulations)
rm -f error.dat errorDominance.dat resultsA_No.res resultsB_No.res ALLresults.res resultsA_HetCat.res resultsA_HetCat1.res resultsA_HetCat2.res resultsB_HetCat.res resultsB_HetCat1.res resultsB_HetCat2.res results.res
for ((j=0 ; j <= 10 ; j++ ))
do
	rm -f resultsA_Rev${ratio[$j]}.res
	rm -f resultsA_HetCat_Rev${ratio[$j]}.res
	rm -f resultsA_HetCat1_Rev${ratio[$j]}.res
	rm -f resultsA_HetCat2_Rev${ratio[$j]}.res
	rm -f resultsB_Rev${ratio[$j]}.res
	rm -f resultsB_HetCat_Rev${ratio[$j]}.res
	rm -f resultsB_HetCat1_Rev${ratio[$j]}.res
	rm -f resultsB_HetCat2_Rev${ratio[$j]}.res
done

ratio=(10)
for ((j=1 ; j <= $nA ; j++ ))
do
		for ((jj=$(($nA+1)) ; jj <= $(($nA+$nB)) ; jj++ ))
		do
					rm -f resultsA_JoinedA${j}B${jj}.res
					rm -f resultsA_JoinedA${j}B${jj}_Rev${ratio[0]}.res
					rm -f resultsB_JoinedA${j}B${jj}.res
					rm -f resultsB_JoinedA${j}B${jj}_Rev${ratio[0]}.res
		done
done


## CREATE EMPTY OUTPUT FILES
touch error.dat errorDominance.dat resultsA_No.res resultsB_No.res ALLresults.res resultsA_HetCat.res resultsA_HetCat1.res resultsA_HetCat2.res resultsB_HetCat.res resultsB_HetCat1.res resultsB_HetCat2.res results.res
ratio=(10 9 8 7 6 5 4 3 2 1 0.5)
for ((j=0 ; j <= 10 ; j++ ))
do
	touch resultsA_Rev${ratio[$j]}.res
	touch resultsA_HetCat_Rev${ratio[$j]}.res
	touch resultsA_HetCat1_Rev${ratio[$j]}.res
	touch resultsA_HetCat2_Rev${ratio[$j]}.res
	touch resultsB_Rev${ratio[$j]}.res
	touch resultsB_HetCat_Rev${ratio[$j]}.res
	touch resultsB_HetCat1_Rev${ratio[$j]}.res
	touch resultsB_HetCat2_Rev${ratio[$j]}.res
done

ratio=(10)
for ((j=1 ; j <= $nA ; j++ ))
do
		for ((jj=$(($nA+1)) ; jj <= $(($nA+$nB)) ; jj++ ))
		do
					touch resultsA_JoinedA${j}B${jj}.res
					touch resultsA_JoinedA${j}B${jj}_Rev${ratio[0]}.res
					touch resultsB_JoinedA${j}B${jj}.res
					touch resultsB_JoinedA${j}B${jj}_Rev${ratio[0]}.res
		done
done


## START ITERATION
for(( NNN = 0 ; NNN < 3334 ; NNN++ ))
do
## RUN1 -- NO REGULATION TO DERMINE THE DOMANAT CYCLE
		kSetGeneration=1
		DominantCycle=-1
		reverse=0
		ratio=(-1)
		HetCat=0
		PozRev=-1
		TModel=1
		JoinedCycles=0
		JoinedPointA=0
		JoinedPointB=0	
		
		gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
		./parameter_generalo_v3  ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[0]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB} #2>/dev/null
		cp k_set k_set_temp
		cp proba_parameters.txt parameters.txt
		
		gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
		./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
		
		gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
		./compactAutCat ${nA} ${nB} #2>/dev/null
		f=$?
		if [ $NNN -eq 0 ];
		then
				cp compactAutCat.c compactAutCat_No.c
				cp proba_parameters.txt parameters_No.txt
				cp res.dat res_No.dat
				dot -Tpng GraphRepresentation.dot > GraphRepresentation_No.png
		fi
		
		printf "$f\n" > TypeOfRunning.res
		DC=`cat TypeOfRunning.res`
		ks=`cat k_set`

		printf "$ks $DC\n" >> results.res
		echo "RUN 1 OK ($NNN)"

## IF DOMINANT CYCLE IS A
		if  [ $DC -eq 1 ];
		then
				echo "Cycle A won!"
				printf "$ks $DC\n" >> resultsA_No.res

## RUN 2 -- CROSS-CATALYSIS TYPE 1 (CC1): E+S = ES, ES+X= ESX; ESX->E+P
				kSetGeneration=2
				HetCat=1
				DominantCycle=$DC
				gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
				./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[0]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
				cp proba_parameters.txt parameters.txt
				cp k_SET k_set_temp				

				gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
				./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
						
				gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
				./compactAutCat ${nA} ${nB} #2>/dev/null
				f=$?
## DRAWING THE GRAPH OF REACTION NETWORK
				if [ $NNN -eq 0 ];
				then
						cp compactAutCat.c compactAutCat_HetCat.c
						cp proba_parameters.txt parameters_HetCat.txt
						cp res.dat res_HetCat.dat
						dot -Tpng GraphRepresentation.dot > GraphRepresentation_HetCat.png
				fi

				printf "$f\n" > TypeOfRunning.res
				TOR=`cat TypeOfRunning.res`
				ks=`cat k_SET`
## SELECT THE RESULTS OF CC FROM CC1 BASED ON k_on2
				kon2P=$((nA+nB+2))
				awk -v kon=$kon2P '{print $kon}' k_SET > kon2.temp
				kon2V=`cat kon2.temp`
				if  [ `echo "$kon2V > 0" | bc` -eq 1 ];
				then
						printf "$ks $TOR\n" >> resultsA_HetCat1.res
				else
						printf "$ks $TOR\n" >> resultsA_HetCat.res
				fi

## IF THE FIRST REACTION WAS CC THEN THE ORDER OF SUBTRAT BINDING WILL CHANGE -> CC2
				if  [ `echo "$kon2V > 0" | bc` -eq 1 ];
				then
## RUN 3 -- CROSS-CATALYSIS TYPE 2: E+X = EX, EX+S= EXS; EXS->E+P
						kSetGeneration=0
						DominantCycle=$DC
						HetCat=2
						gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[0]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
						cp proba_parameters.txt parameters.txt
						cp k_SET2 k_set_temp						

						gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
						
						gcc -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
						./compactAutCat ${nA} ${nB} #2>/dev/null
						f=$?
						if [ $NNN -eq 0 ];
						then
								cp compactAutCat.c compactAutCat_HetCat2.c
								cp proba_parameters.txt parameters_HetCat2.txt
								cp res.dat res_HetCat2.dat
								dot -Tpng GraphRepresentation.dot > GraphRepresentation_HetCat2.png
						fi

						printf "$f\n" > TypeOfRunning.res
						TOR=`cat TypeOfRunning.res`
						ks=`cat k_SET2`

						printf "$ks $TOR\n" >> resultsA_HetCat2.res
				fi
				ratio=(10 9 8 7 6 5 4 3 2 1 0.5)
## EFFECT OF REVERSE REACTION AT THE POZREV POZITION WHITOUTH CROSS-CATALYSIS AND WITH CC, CC1 AND CC2
				for ((j=0 ; j <= 10 ; j++ ))
				do
## RUN 4 -- REVERSE ON LAST REACTION WITHOUT CROSS-CATALYSIS
						kSetGeneration=0
						DominantCycle=$DC
						HetCat=0
						reverse=1
						PozRev=$(( ${nA}))
						gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[$j]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
						cp proba_parameters.txt parameters.txt
						cp k_SET2 k_set_temp
      
						gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
						
						gcc -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
						./compactAutCat ${nA} ${nB} #2>/dev/null
						f=$?
						if [ $NNN -eq 0 ];
						then
								cp compactAutCat.c compactAutCat_Rev${ratio[$j]}.c
								cp proba_parameters.txt parameters_Rev${ratio[$j]}.txt
								cp res.dat res_Rev${ratio[$j]}.dat
								dot -Tpng GraphRepresentation.dot > GraphRepresentation_Rev${ratio[$j]}.png
						fi
						
						printf "$f\n" > TypeOfRunning.res
						TOR=`cat TypeOfRunning.res`
						ks=`cat k_SET2`

						printf "$ks $TOR\n" >> resultsA_Rev${ratio[$j]}.res

## RUN 5 -- REVERSE ON LAST REACTION AND CROSS-CATALYSIS TYPE 1: E+S = ES, ES+X= ESX; ESX->E+P
						kSetGeneration=0
						DominantCycle=$DC
						reverse=1
						HetCat=1
						PozRev=$(( ${nA}))
						gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[$j]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
						cp proba_parameters.txt parameters.txt
						cp k_SET2 k_set_temp						

						gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
						
						gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
						./compactAutCat ${nA} ${nB} #2>/dev/null
						f=$?
						if [ $NNN -eq 0 ];
						then
								cp compactAutCat.c compactAutCat_HetCat_Rev${ratio[$j]}.c
								cp proba_parameters.txt parameters_HetCat_Rev${ratio[$j]}.txt
								cp res.dat res_HetCat_Rev${ratio[$j]}.dat
								dot -Tpng GraphRepresentation.dot > GraphRepresentation_HetCat_Rev${ratio[$j]}.png
						fi

						printf "$f\n" > TypeOfRunning.res
						TOR=`cat TypeOfRunning.res`
						ks=`cat k_SET2`

						kon2P=$((nA+nB+2))
						awk -v kon=$kon2P '{print $kon}' k_SET2 > kon2.temp
						kon2V=`cat kon2.temp`

			
						if  [ `echo "$kon2V > 0" | bc` -eq 1 ];
						then
								printf "$ks $TOR\n" >> resultsA_HetCat1_Rev${ratio[$j]}.res
						else
								printf "$ks $TOR\n" >> resultsA_HetCat_Rev${ratio[$j]}.res
						fi
			
						if  [ `echo "$kon2V > 0" | bc` -eq 1 ];
						then
## RUN 6 -- REVERSE ON LAST REACTION AND CROSS-CATALYSIS TYPE 2: E+X = EX, EX+S= EXS; EXS->E+P
								kSetGeneration=0
								DominantCycle=$DC
								HetCat=2
								reverse=1
								PozRev=$(( ${nA}))
								gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[$j]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
								cp proba_parameters.txt parameters.txt
								cp k_SET2 k_set_temp

								gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
								
								gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
								./compactAutCat ${nA} ${nB} #2>/dev/null
								f=$?
								if [ $NNN -eq 0 ];
								then
										cp compactAutCat.c compactAutCat_HetCat2_Rev${ratio[$j]}.c
										cp proba_parameters.txt parameters_HetCat2_Rev${ratio[$j]}.txt
										cp res.dat res_HetCat2_Rev${ratio[$j]}.dat
										dot -Tpng GraphRepresentation.dot > GraphRepresentation_HetCat2_Rev${ratio[$j]}.png
								fi								
								

								printf "$f\n" > TypeOfRunning.res
								TOR=`cat TypeOfRunning.res`
								ks=`cat k_SET2`

								printf "$ks $TOR\n" >> resultsA_HetCat2_Rev${ratio[$j]}.res

						fi
				done
## RUN 7-15 -- JOINED Cycle A and Cycle B at Ai=Bj
				for ((j=1 ; j <= $nA ; j++ ))
				do
						for ((jj=$(($nA+1)) ; jj <= $(($nA+$nB)) ; jj++ ))
						do
								kSetGeneration=0
								DominantCycle=$DC
								HetCat=0
								reverse=0
								PozRev=-1 #$(( ${nA}-1 ))
								JoinedCycles=1
								JoinedPointA=$j
								JoinedPointB=$jj
								ratio=(-1)
								gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[0]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
								cp proba_parameters.txt parameters.txt
								cp k_set k_set_temp
								
								gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
								
								gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
								./compactAutCat ${nA} ${nB} #2>/dev/null
								f=$?
								if [ $NNN -eq 0 ];
								then
										cp compactAutCat.c compactAutCat_JoinedA${j}B${jj}.c
										cp proba_parameters.txt parameters_JoinedA${j}B${jj}.txt
										cp res.dat res_JoinedA${j}B${jj}.dat
										dot -Tpng GraphRepresentation.dot > GraphRepresentation_JoinedA${j}B${jj}.png
								fi								
								
								printf "$f\n" > TypeOfRunning.res
								TOR=`cat TypeOfRunning.res`
								ks=`cat k_set`
								printf "$ks $TOR\n" >> resultsA_JoinedA${j}B${jj}.res
						done
				done
				JoinedCycles=0
				JoinedPointA=0
				JoinedPointB=0
## RUN 16-24 -- JOINED Cycle A and Cycle B at Ai=Bj AND REVERSE
				for ((j=1 ; j <= $nA ; j++ ))
				do
						for ((jj=$(($nA+1)) ; jj <= $(($nA+$nB)) ; jj++ ))
						do
								kSetGeneration=0
								DominantCycle=$DC
								HetCat=0
								reverse=1
								PozRev=$(( ${nA}))
								JoinedCycles=1
								JoinedPointA=$j
								JoinedPointB=$jj
								ratio=(10)
								gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[0]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
								cp proba_parameters.txt parameters.txt
								cp k_set k_set_temp
								
								gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
								
								gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
								./compactAutCat ${nA} ${nB} #2>/dev/null
								f=$?
								if [ $NNN -eq 0 ];
								then
										cp compactAutCat.c compactAutCat_JoinedA${j}B${jj}_Rev${ratio[0]}.c
										cp proba_parameters.txt parameters_JoinedA${j}B${jj}_Rev${ratio[0]}.txt
										cp res.dat res_JoinedA${j}B${jj}_Rev${ratio[0]}.dat
										dot -Tpng GraphRepresentation.dot > GraphRepresentation_JoinedA${j}B${jj}_Rev${ratio[0]}.png
								fi								
								
								printf "$f\n" > TypeOfRunning.res
								TOR=`cat TypeOfRunning.res`
								ks=`cat k_set`
								printf "$ks $TOR\n" >> resultsA_JoinedA${j}B${jj}_Rev${ratio[0]}.res
						done
				done
				JoinedCycles=0
				JoinedPointA=0
				JoinedPointB=0
## IF THE DOMINANT CYCLE IS B
		elif [ $DC -eq 2 ];
		then
				echo "Cycle B won!"
				printf "$ks $DC\n" >> resultsB_No.res
				
## RUN 2 -- CROSS-CATALYSIS TYPE 1: E+S = ES, ES+X= ESX; ESX->E+P
				kSetGeneration=2
				HetCat=1
				DominantCycle=$DC
				ratio=(-1)
				gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
				./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[0]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
				cp proba_parameters.txt parameters.txt
				cp k_SET k_set_temp				

				gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
				./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
						
				gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
				./compactAutCat ${nA} ${nB} #2>/dev/null
				f=$?
				if [ $NNN -eq 0 ];
				then
						cp compactAutCat.c compactAutCat_HetCat.c
						cp proba_parameters.txt parameters_HetCat.txt
						cp res.dat res_HetCat.dat
						dot -Tpng GraphRepresentation.dot > GraphRepresentation_HetCat.png
				fi
				
				printf "$f\n" > TypeOfRunning.res
				TOR=`cat TypeOfRunning.res`
				ks=`cat k_SET`

				kon2P=$((nA+nB+2))
				awk -v kon=$kon2P '{print $kon}' k_SET > kon2.temp
				kon2V=`cat kon2.temp`
				if  [ `echo "$kon2V > 0" | bc` -eq 1 ];
				then
						printf "$ks $TOR\n" >> resultsB_HetCat1.res
				else
						printf "$ks $TOR\n" >> resultsB_HetCat.res
				fi

				if  [ `echo "$kon2V > 0" | bc` -eq 1 ];
				then
## RUN 3 -- CROSS-CATALYSIS TYPE 2: E+X = EX, EX+S= EXS; EXS->E+P
						kSetGeneration=0
						DominantCycle=$DC
						HetCat=2
						gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[0]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
						cp proba_parameters.txt parameters.txt
						cp k_SET2 k_set_temp						

						gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
						
						gcc -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
						./compactAutCat ${nA} ${nB} #2>/dev/null
						f=$?
						if [ $NNN -eq 0 ];
						then
								cp compactAutCat.c compactAutCat_HetCat2.c
								cp proba_parameters.txt parameters_HetCat2.txt
								cp res.dat res_HetCat2.dat
								dot -Tpng GraphRepresentation.dot > GraphRepresentation_HetCat2.png
						fi
						
						printf "$f\n" > TypeOfRunning.res
						TOR=`cat TypeOfRunning.res`
						ks=`cat k_SET2`

						printf "$ks $TOR\n" >> resultsB_HetCat2.res

				fi
				ratio=(10 9 8 7 6 5 4 3 2 1 0.5)
				#ratio=(10)
				for ((j=0 ; j <= 10 ; j++ ))
				do
## RUN 4 -- REVERSE ON LAST REACTION
						kSetGeneration=0
						DominantCycle=$DC
						HetCat=0
						reverse=1
						PozRev=$(( ${nB}))
						gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[$j]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
						cp proba_parameters.txt parameters.txt
						cp k_SET2 k_set_temp
      
						gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
						
						gcc -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
						./compactAutCat ${nA} ${nB} #2>/dev/null
						f=$?
						if [ $NNN -eq 0 ];
						then
								cp compactAutCat.c compactAutCat_Rev${ratio[$j]}.c
								cp proba_parameters.txt parameters_Rev${ratio[$j]}.txt
								cp res.dat res_Rev${ratio[$j]}.dat
								dot -Tpng GraphRepresentation.dot > GraphRepresentation_Rev${ratio[$j]}.png
						fi

						printf "$f\n" > TypeOfRunning.res
						TOR=`cat TypeOfRunning.res`
						ks=`cat k_SET2`

						printf "$ks $TOR\n" >> resultsB_Rev${ratio[$j]}.res

## RUN 5 -- REVERSE ON LAST REACTION AND CROSS-CATALYSIS TYPE 1: E+S = ES, ES+X= ESX; ESX->E+P
						kSetGeneration=0
						DominantCycle=$DC
						reverse=1
						HetCat=1
						PozRev=$(( ${nB}))
						gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[$j]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
						cp proba_parameters.txt parameters.txt
						cp k_SET2 k_set_temp						

						gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
						./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
						
						gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
						./compactAutCat ${nA} ${nB} #2>/dev/null
						f=$?
						if [ $NNN -eq 0 ];
						then
								cp compactAutCat.c compactAutCat_HetCat_Rev${ratio[$j]}.c
								cp proba_parameters.txt parameters_HetCat_Rev${ratio[$j]}.txt
								cp res.dat res_HetCat_Rev${ratio[$j]}.dat
								dot -Tpng GraphRepresentation.dot > GraphRepresentation_HetCat_Rev${ratio[$j]}.png
						fi
						

						printf "$f\n" > TypeOfRunning.res
						TOR=`cat TypeOfRunning.res`
						ks=`cat k_SET2`

						kon2P=$((nA+nB+2))
						awk -v kon=$kon2P '{print $kon}' k_SET2 > kon2.temp
						kon2V=`cat kon2.temp`
			
						if  [ `echo "$kon2V > 0" | bc` -eq 1 ];
						then
								printf "$ks $TOR\n" >> resultsB_HetCat1_Rev${ratio[$j]}.res
						else
								printf "$ks $TOR\n" >> resultsB_HetCat_Rev${ratio[$j]}.res
						fi
			
						if  [ `echo "$kon2V > 0" | bc` -eq 1 ];
						then
## RUN 6 -- REVERSE ON LAST REACTION AND CROSS-CATALYSIS TYPE 2: E+X = EX, EX+S= EXS; EXS->E+P
								kSetGeneration=0
								DominantCycle=$DC
								HetCat=2
								reverse=1
								PozRev=$(( ${nB}))
								gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[$j]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
								cp proba_parameters.txt parameters.txt
								cp k_SET2 k_set_temp

								gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
								
								gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
								./compactAutCat ${nA} ${nB} #2>/dev/null
								f=$?
								if [ $NNN -eq 0 ];
								then
										cp compactAutCat.c compactAutCat_HetCat2_Rev${ratio[$j]}.c
										cp proba_parameters.txt parameters_HetCat2_Rev${ratio[$j]}.txt
										cp res.dat res_HetCat2_Rev${ratio[$j]}.dat
										dot -Tpng GraphRepresentation.dot > GraphRepresentation_HetCat2_Rev${ratio[$j]}.png
								fi								
								
								printf "$f\n" > TypeOfRunning.res
								TOR=`cat TypeOfRunning.res`
								ks=`cat k_SET2`
								printf "$ks $TOR\n" >> resultsB_HetCat2_Rev${ratio[$j]}.res
						fi
				done
## RUN 7-15 -- JOINED Cycle A and Cycle B at Ai=Bj
				for ((j=1 ; j <= $nA ; j++ ))
				do
						for ((jj=$(($nA+1)) ; jj <= $(($nA+$nB)) ; jj++ ))
						do
								kSetGeneration=0
								DominantCycle=$DC
								HetCat=0
								reverse=0
								PozRev=-1 #$(( ${nA}-1 ))
								JoinedCycles=1
								JoinedPointA=$j
								JoinedPointB=$jj
								ratio=(-1)
								gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[0]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
								cp proba_parameters.txt parameters.txt
								cp k_set k_set_temp
								
								gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
								
								gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
								./compactAutCat ${nA} ${nB} #2>/dev/null
								f=$?
								if [ $NNN -eq 0 ];
								then
										cp compactAutCat.c compactAutCat_JoinedA${j}B${jj}.c
										cp proba_parameters.txt parameters_JoinedA${j}B${jj}.txt
										cp res.dat res_JoinedA${j}B${jj}.dat
										dot -Tpng GraphRepresentation.dot > GraphRepresentation_JoinedA${j}B${jj}.png
								fi								
								
								printf "$f\n" > TypeOfRunning.res
								TOR=`cat TypeOfRunning.res`
								ks=`cat k_set`
								printf "$ks $TOR\n" >> resultsB_JoinedA${j}B${jj}.res
						done
				done
				JoinedCycles=0
				JoinedPointA=0
				JoinedPointB=0
## RUN 16-24 -- JOINED Cycle A and Cycle B at Ai=Bj AND REVERSE
				for ((j=1 ; j <= $nA ; j++ ))
				do
						for ((jj=$(($nA+1)) ; jj <= $(($nA+$nB)) ; jj++ ))
						do
								kSetGeneration=0
								DominantCycle=$DC
								HetCat=0
								reverse=1
								PozRev=$(( ${nB}))
								JoinedCycles=1
								JoinedPointA=$j
								JoinedPointB=$jj
								ratio=(10)
								gcc -Wall -pedantic -Ofast -o parameter_generalo_v3 parameter_generalo_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parameter_generalo_v3 ${nA} ${nB} ${kSetGeneration} ${DominantCycle} ${InitRandNum} ${reverse} ${ratio[0]} ${HetCat} ${PozRev} ${TModel} ${JoinedCycles} ${JoinedPointA} ${JoinedPointB}
								cp proba_parameters.txt parameters.txt
								cp k_set k_set_temp
								
								gcc -pedantic -Ofast -o parsing_v3 parsing_v3.c -L/usr/local/lib/ -lgsl -lgslcblas -lm 2>/dev/null
								./parsing_v3 parameters.txt ${nA} ${nB} #2>/dev/null
								
								gcc -Wall -pedantic -Ofast -o compactAutCat compactAutCat.c -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -I/usr/local/include -lm 2>/dev/null
								./compactAutCat ${nA} ${nB} #2>/dev/null
								f=$?
								if [ $NNN -eq 0 ];
								then
										cp compactAutCat.c compactAutCat_JoinedA${j}B${jj}_Rev${ratio[0]}.c
										cp proba_parameters.txt parameters_JoinedA${j}B${jj}_Rev${ratio[0]}.txt
										cp res.dat res_JoinedA${j}B${jj}_Rev${ratio[0]}.dat
										dot -Tpng GraphRepresentation.dot > GraphRepresentation_JoinedA${j}B${jj}_Rev${ratio[0]}.png
								fi								
								
								printf "$f\n" > TypeOfRunning.res
								TOR=`cat TypeOfRunning.res`
								ks=`cat k_set`

								printf "$ks $TOR\n" >> resultsB_JoinedA${j}B${jj}_Rev${ratio[0]}.res
						done
				done
				JoinedCycles=0
				JoinedPointA=0
				JoinedPointB=0

		else
				echo "Nobody won!"
				printf "$ks $DC\n" >> errorDominance.dat
		fi
done

		
