```bash
#!/bin/bash


filename="$(echo $1 | tr -d '.c')"
#echo $filename

gcc Param.c # to create the file
./a.out | sort -n | uniq > Parameter1.txt # to create the file 




for par in $(cat Parameter1.txt)
do
#mkdir ${par}_${filename}_dir 
 #iter=PARA

 cat $1 | sed s/-9084099491-PARA/-9084099491-${par}/ > Ran_id${par}/${filename}_${par}.c # seed for random number
 

  
 gcc Ran_id${par}/${filename}_${par}.c -lm -O2 -o Ran_id${par}/${filename}_${par}.out
 #./${filename}_${par}.out
 cd Ran_id${par}/
 nohup ./${filename}_${par}.out &> ${filename}_${par}.ave &
 
 cd ..
 
 
done
```
