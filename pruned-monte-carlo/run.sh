#!/bin/bash
#LOOP=10
GRAPH=("soc-Epinions1.txt" "soc-LiveJournal1.txt" "dblp.txt")
GRAPH_LENGTH=${#GRAPH[@]}

KSEED=("10" "20" "30" "40" "50")
KSEED_LENGTH=${#KSEED[@]}

PRO=("0.01" "0.025" "0.05" "0.075" "0.1")
PRO_LENGHT=${#PRO[@]}

#Graph
for (( i=0; i<${GRAPH_LENGTH}; i++ ));
do
        echo -e "Performing ${GRAPH[$i]}"
        #Seed size
        for (( j=0; j<${KSEED_LENGTH}; j++ ));
        do
                echo -e "\tk=${KSEED[$j]}---"
                #Propagation
                for (( k=0; k<${PRO_LENGHT}; k++ ));
                do
                        echo -e "\t\tP=${PRO[$k]}---"
                        #R  unnecessary
                        ./benchmark_estimate ../data/${GRAPH[$i]} ${KSEED[$j]} 10 ${PRO[$k]}
                        echo
                done
                echo
        done
        echo
done
