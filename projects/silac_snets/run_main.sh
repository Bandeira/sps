#!/bin/bash                                                                                                                                    
#                                                                                                                                              
#$ -cwd                                                                                                                                        
#$ -j y                                                                                                                                        
#$ -S /bin/bash                                                                                                                                
#                                                                                                                                              
export LD_LIBRARY_PATH=../../install/linux-g++-64/lib

echo "../../trunk/ExecFramework/main_specnets" $1 $2 $3 $4 $5 $6 $7 $8 $9

../../trunk/ExecFramework/main_specnets $1 $2 $3 $4 $5 $6 $7 $8 $9
