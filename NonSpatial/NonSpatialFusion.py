import numpy as np
import math
from numba import njit
from numba import jit
import time
import sys
import random
import os



#Functions
@jit
def countPopulation(ListePopulation):
    countPop=0
    for i in range (0,len(ListePopulation)):
        countPop=countPop+ListePopulation[i][0]
    return countPop



@jit
def countSpecies(ListePopulation):
    countS=0
    for i in range (0,len(ListePopulation)):
        if (ListePopulation[i][0]>0):
            countS=countS+1
    return countS



@njit
def convertBinaire(x):
    n=len(x)
    res=0
    for i in range (0,n):
        res=res+x[i]*(2**i)
    return res


@jit
def ComputeShanon(ListePopulation):
    res=0
    Tot=0
    prop=0
    for i in range (0,len(ListePopulation)):
        Tot=Tot+ListePopulation[i][0]

    for i in range (0,len(ListePopulation)):
        if (ListePopulation[i][0]>0):
            prop=float(ListePopulation[i][0])/Tot
            res=res-prop*np.log(prop)

    return res


@jit
def ComputeIndex(ListePopulation,q):
    res=0
    if (q==1):
        res=ComputeShanon(ListePopulation)
    else:

        Tot=0
        prop=0;
        for i in range (0,len(ListePopulation)):
            Tot=Tot+ListePopulation[i][0]

        for i in range (0,len(ListePopulation)):
            prop=float(ListePopulation[i][0])/Tot
            res=res+prop**q
        res=res**(float(1)/(1-q))

    return res




@jit
def Mutation(x):
    res=np.zeros(len(x))
    for i in range (0,len(x)):
        res[i]=x[i]
    PartieNulle=np.zeros(len(x),dtype=int)
    count=0
    for i in range(0,len(res)):
        if (res[i]==0):
            PartieNulle[count]=i
            count=count+1
    if (count>0):
        m=np.random.randint(0,count-1)
        res[PartieNulle[m]]=1
    return res




@jit
def chooseElement(ListePopulation):
    countPop=0
    res=0
    for i in range (0,len(ListePopulation)):
        countPop=countPop+ListePopulation[i][0]

    if (countPop==1):
        y=1
    else:
        y=np.random.randint(1,countPop+1)
    countPop2=0
    count=0
    found=0
    while (count<len(ListePopulation) and found==0):
        if (ListePopulation[count][0]>0):
            countPop2=countPop2+ListePopulation[count][0]
        if (y<=countPop2 and ListePopulation[count][0]>0):
            found=1
            res=count
        count=count+1
    if (found==0):
        print("Element not found!!!!!")
    return res


@jit
def cleanData(ListePopulation):
    j=0
    while(j<len(ListePopulation)):
        if (ListePopulation[j][0]==0):
            del(ListePopulation[j])
        else:
            ListePopulation[j][2]=0
            ListePopulation[j][5]=0
            j=j+1


@njit
def minimum(x1,x2):
    res=x1
    if (x2<x1):
        res=x2
    return res




@njit
def fusionVect(G1,G2):
    Off1=np.zeros(len(G1))
    Off2=np.zeros(len(G1))
    
    for i in range (0,len(G1)):
        y=np.random.randint(0,2)
        if (y==0):
            Off1[i]=G1[i]
            Off2[i]=G2[i]
        else:
            Off1[i]=G2[i]
            Off2[i]=G1[i]


    return (Off1,Off2)


@njit 
def numbaSum(Genotype):

    res=np.sum(Genotype)

    return res



@jit
def MaxScore(ListePopulation):
    res=0
    for i in range (0,len(ListePopulation)):
        scorei=numbaSum(ListePopulation[i][1])
        if (res<scorei):
            res=scorei
    return res


start = time.time()


#Number of considered genes
# Ngenes=300 #Number of driver genes
Ngenes=50 #Number of driver genes
Genotype=np.zeros(Ngenes)



#Parameters entered as a terminal input
KCStr=sys.argv[2]
KC=int(KCStr)
NgenerationsMaxStr=sys.argv[3]
NgenerationsMax=int(NgenerationsMaxStr)
pmStr=sys.argv[4]
pm=float(pmStr)
pfStr=sys.argv[5]
pf=float(pfStr)
growthRateStr=sys.argv[6]
growthRate=float(growthRateStr)
deathRate=growthRate
DT=1
# DT=0.1



GenotypeTemporaire=np.zeros(Ngenes)

Number=np.zeros(NgenerationsMax) #Number of genotypes per time step
Shanon=np.zeros(NgenerationsMax) #Shanon index per time step
Simpson=np.zeros(NgenerationsMax) #Simpson index per time step
Score=np.zeros(NgenerationsMax) #Maximal amount of mutations per time step
TotalCells=np.zeros(NgenerationsMax) #Total number of cells per time step
FinalList=[]


# Results folder path
directory_path='Results/CC'+str(KC)+'/Neutral/LogisticFusion/g'+str(growthRate)+'/mu'+str(pm)+'/pf'+str(pf)+'/'
directory=os.path.dirname(directory_path)


if not os.path.exists(directory):
    os.makedirs(directory)


#Main Function
@jit
def ModelRun(NgenerationsMax,DT,s,KC,pm):

    np.random.seed(int(s))
    ExtinctSexSpecie=0
    OffspringM=0
    ListePop=[]
    ListeExtincted=[]
    genotypesCounts=0
    ListePop.append([1,Genotype,0,0,genotypesCounts,0,genotypesCounts])


    #Creating and opening results file
    FileName1=directory_path+'Lineage'+str(s)+'.txt'
    FileName2=directory_path+'Total'+str(s)+'.txt'
    FileName3=directory_path+'NumberMut'+str(s)+'.txt'
    FileName4=directory_path+'Shanon'+str(s)+'.txt'
    FileName5=directory_path+'Simpson'+str(s)+'.txt'
    FileName6=directory_path+'Score'+str(s)+'.txt'


    f1=open(FileName1,'w')
    f1.write("Size Id Time Ancestor\n")
    f2=open(FileName2,'w')
    f3=open(FileName3,'w')
    f4=open(FileName4,'w')
    f5=open(FileName5,'w')
    f6=open(FileName6,'w')


    #Time loop
    for l in range (0,NgenerationsMax):
        print('Generation=',l)
        newD=0
        TotalPopulation=countPopulation(ListePop)

        #Population loop
        for j in range(0,len(ListePop)):
            nombreRepresentants=ListePop[j][0]
            if (nombreRepresentants>0):
                newCells=np.random.poisson(growthRate*nombreRepresentants*DT) #Newborns
                newM=np.random.binomial(newCells,Ngenes*pm) #New mutants
                if (newM>newCells):
                    newM=newCells
                
                newD=np.random.poisson(deathRate*nombreRepresentants*TotalPopulation*DT/KC) #New dead cells

                countDeads=0
                newD1=0
                newD2=0
                
                if (newD<nombreRepresentants+newCells):
                    while (countDeads<newD):
                        m=np.random.randint(1,nombreRepresentants+newCells)
                        if (m<nombreRepresentants+newCells-newM and newD1<nombreRepresentants+newCells-newM):
                            newD1=newD1+1
                        else:
                            newD2=newD2+1
                        countDeads=countDeads+1
                    if (ListePop[j][0]+newCells-newM-newD1>0):
                        ListePop[j][0]=ListePop[j][0]+newCells-newM-newD1
                    else:
                        ListePop[j][0]=0
                    if (newM-newD2>0):
                        ListePop[j][2]=newM-newD2
                    else:
                        ListePop[j][2]=0
                else:
                    ListePop[j][0]=0
                    ListePop[j][2]=0

                if (ListePop[j][0]==0):
                    NewElement=[1,ListePop[j][1],ListePop[j][2],l*DT,ListePop[j][4],ListePop[j][5],ListePop[j][6]]
                    ListeExtincted.append(NewElement)


        j=0
        Ncurrent=len(ListePop)
        while(j<Ncurrent):
            for k in range(0,ListePop[j][2]):
                GenotypeTemporaire=Mutation(ListePop[j][1])
                exist=0
                count=0
                while (count<len(ListePop) and exist==0):
                        
                    if (convertBinaire(ListePop[count][1])-convertBinaire(GenotypeTemporaire)==0):
                        ListePop[count][0]=ListePop[count][0]+1
                        exist=1
                    count=count+1
                if (exist==0):
                    genotypesCounts=genotypesCounts+1
                    NewElement=[1,GenotypeTemporaire,0,l*DT,ListePop[j][6],0,genotypesCounts]
                    ListePop.append(NewElement)
            ListePop[j][2]=0
            j=j+1
        
        #Hybrids formation
        TotalPopulation=countPopulation(ListePop)
        newH=np.random.poisson(pf*TotalPopulation*DT) #New hybrids
        newH=minimum(newH,int(TotalPopulation/2))


        for j in range (0,newH):
            y=np.random.randint(1,TotalPopulation+1)
            countPop2=0
            count=0
            found=0
            
            while (count<len(ListePop) and found==0):
                if (ListePop[count][0]>0):
                    countPop2=countPop2+ListePop[count][0]
                    if (y<=countPop2):
                        found=1
                        ListePop[count][5]=ListePop[count][5]+1
                        ListePop[count][0]=ListePop[count][0]-1
                        TotalPopulation=TotalPopulation-1

                    if (ListePop[count][0]==0):
                        NewElement=[1,ListePop[count][1],ListePop[count][2],l*DT,ListePop[count][4],ListePop[count][5],ListePop[count][6]]
                        ListeExtincted.append(NewElement)
                count=count+1


        Ncurrent=len(ListePop)
        TotalPopulation=countPopulation(ListePop)

        if (TotalPopulation>0):
            for j in range(0,Ncurrent):
                for k in range(0,ListePop[j][5]):
                    Genotype1=ListePop[j][1]

                    neighbor=chooseElement(ListePop)
                    ListePop[neighbor][0]=ListePop[neighbor][0]-1
                    if (ListePop[neighbor][0]==0):
                        NewElement=[1,ListePop[neighbor][1],ListePop[neighbor][2],l*DT,ListePop[neighbor][4],ListePop[neighbor][5],ListePop[neighbor][6]]
                        ListeExtincted.append(NewElement)

                    Genotype2=ListePop[neighbor][1]

                        
                    (Off1,Off2)=fusionVect(Genotype1,Genotype2) #Recombination

                    exist=0
                    count=0
                    while (count<len(ListePop) and exist==0):
                        if (convertBinaire(ListePop[count][1])-convertBinaire(Off1)==0):

                            ListePop[count][0]=ListePop[count][0]+1
                            exist=1
                        count=count+1
                    if (exist==0):
                        genotypesCounts=genotypesCounts+1


                        NewElement=[1,Off1,0,l*DT,ListePop[j][6],0,genotypesCounts] #New genotype!
                        ListePop.append(NewElement)

                    
                    exist=0
                    count=0
                    while (count<len(ListePop) and exist==0):
                        if (convertBinaire(ListePop[count][1])-convertBinaire(Off2)==0):

                            ListePop[count][0]=ListePop[count][0]+1
                            exist=1
                        count=count+1
                    if (exist==0):
                        genotypesCounts=genotypesCounts+1

                        NewElement=[1,Off2,0,l*DT,ListePop[j][6],0,genotypesCounts] #New genotype!
                        ListePop.append(NewElement)

                    
                ListePop[j][5]=0

        else:
            for j in range(0,Ncurrent):
                for k in range(0,ListePop[j][5]):
                    ListePop[j][0]=ListePop[j][0]+1


        cleanData(ListePop)

        

        Number[l]=countSpecies(ListePop)
        TotalCells[l]=countPopulation(ListePop)
        Shanon[l]=ComputeShanon(ListePop)
        Simpson[l]=ComputeIndex(ListePop,2)
        Score[l]=MaxScore(ListePop)

        

        #Writting results
        if (l%50==0 or l==NgenerationsMax-1):


                for i in range(0,len(ListeExtincted)):
                        f1.write("%i %i %5.1f %i\n" % (1,ListeExtincted[i][6],l*DT,ListeExtincted[i][4]))

                for i in range(0,len(ListePop)):
                        f1.write("%i %i %5.1f %i\n" % (ListePop[i][0],ListePop[i][6],l*DT,ListePop[i][4]))

                f2.write("%i %i\n" % (TotalCells[l],l*DT))
                f3.write("%i %i\n" % (Number[l],l*DT))
                f4.write("%i %i\n" % (Shanon[l],l*DT))
                f5.write("%i %i\n" % (Simpson[l],l*DT))
                f6.write("%i %i\n" % (Score[l],l*DT))

                ListeExtincted=[]


    #Closing file    
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    f6.close()
    


s=sys.argv[1]


ModelRun(NgenerationsMax,DT,s,KC,pm)

end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))