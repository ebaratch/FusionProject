import numpy as np
import math
from numba import njit
from numba import jit
import time
import sys
import random
import os



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
def TargetedMutation(x,j):
    res=np.zeros(len(x))
    for i in range (0,len(x)):
        res[i]=x[i]
    PartieNulle=np.zeros(len(x),dtype=int)
    count=0
    if (res[j]==0):
        res[j]=1;

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
    #print("countPop=",countPop)
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

Ngenes=300 #Number of driver genes
# Ngenes=50 #Number of driver genes
Genotype=np.zeros(Ngenes)
GenotypeM1=np.zeros(Ngenes)
GenotypeM2=np.zeros(Ngenes)
GenotypeM3=np.zeros(Ngenes)
GenotypeM4=np.zeros(Ngenes)
GenotypeM5=np.zeros(Ngenes)


    
GenotypeM1[0]=1
GenotypeM2[1]=1
GenotypeM3[2]=1
GenotypeM4[3]=1
GenotypeM5[4]=1


KCStr=sys.argv[2]
KC=int(KCStr)
pmStr=sys.argv[3]
pm=float(pmStr)
pfStr=sys.argv[4]
pf=float(pfStr)
growthRateStr=sys.argv[5]
growthRate=float(growthRateStr)
diversityStr=sys.argv[6]
diversity=int(diversityStr)
deathRate=growthRate
NgenerationsMaxStr=sys.argv[7]
NgenerationsMax=int(NgenerationsMaxStr)
DTStr=sys.argv[8]
DT=float(DTStr)




GenotypeTemporaire=np.zeros(Ngenes)

Number=np.zeros(NgenerationsMax)
Shanon=np.zeros(NgenerationsMax)
Simpson=np.zeros(NgenerationsMax)
Score=np.zeros(NgenerationsMax)
TotalCells=np.zeros(NgenerationsMax)
FinalList=[]


directory_path='Results/CC'+str(KC)+'/Neutral/Diversity/'+str(diversity)+'/LogisticFusion/g'+str(growthRate)+'/mu'+str(pm)+'/pf'+str(pf)+'/'
directory=os.path.dirname(directory_path)


if not os.path.exists(directory):
    os.makedirs(directory)


@jit
def ModelRun(NgenerationsMax,DT,s,KC,pm,pf,diversity):

    np.random.seed(int(s))
    ExtinctSexSpecie=0
    OffspringM=0
    ListePop=[]
    ListeExtincted=[]
    genotypesCounts=0
    ListePop.append([1,Genotype,0,0,genotypesCounts,0,genotypesCounts])
    for i in range (0,diversity):
        ListePop.append([1,TargetedMutation(Genotype,i),0,0,genotypesCounts,0,genotypesCounts])


    #Creating and opening results file

    FileName1Ini=directory_path+'TotalIni'+str(s)+'.txt'
    FileName2Ini=directory_path+'NumberMutIni'+str(s)+'.txt'
    FileName3Ini=directory_path+'ShanonIni'+str(s)+'.txt'
    FileName4Ini=directory_path+'SimpsonIni'+str(s)+'.txt'
    FileName5Ini=directory_path+'ScoreIni'+str(s)+'.txt'


    FileName1First=directory_path+'TotalFirst'+str(s)+'.txt'
    FileName2First=directory_path+'NumberMutFirst'+str(s)+'.txt'
    FileName3First=directory_path+'ShanonFirst'+str(s)+'.txt'
    FileName4First=directory_path+'SimpsonFirst'+str(s)+'.txt'
    FileName5First=directory_path+'ScoreFirst'+str(s)+'.txt'


    FileName1=directory_path+'Total'+str(s)+'.txt'
    FileName2=directory_path+'NumberMut'+str(s)+'.txt'
    FileName3=directory_path+'Shanon'+str(s)+'.txt'
    FileName4=directory_path+'Simpson'+str(s)+'.txt'
    FileName5=directory_path+'Score'+str(s)+'.txt'



    f1Ini=open(FileName1Ini,'w')
    f2Ini=open(FileName2Ini,'w')
    f3Ini=open(FileName3Ini,'w')
    f4Ini=open(FileName4Ini,'w')
    f5Ini=open(FileName5Ini,'w')



    f1First=open(FileName1First,'w')
    f2First=open(FileName2First,'w')
    f3First=open(FileName3First,'w')
    f4First=open(FileName4First,'w')
    f5First=open(FileName5First,'w')



    f1=open(FileName1,'w')
    f2=open(FileName2,'w')
    f3=open(FileName3,'w')
    f4=open(FileName4,'w')
    f5=open(FileName5,'w')


    NumberIni=countSpecies(ListePop)
    TotalCellsIni=countPopulation(ListePop)
    ShanonIni=ComputeShanon(ListePop)
    SimpsonIni=ComputeIndex(ListePop,2)
    ScoreIni=MaxScore(ListePop)



    f1Ini.write("%i\n" % (TotalCellsIni))
    f2Ini.write("%i\n" % (NumberIni))
    f3Ini.write("%i\n" % (ShanonIni))
    f4Ini.write("%i\n" % (SimpsonIni))
    f5Ini.write("%i\n" % (ScoreIni))


    f1Ini.close()
    f2Ini.close()
    f3Ini.close()
    f4Ini.close()
    f5Ini.close()


    for l in range (0,NgenerationsMax):
        print('Generation=',l)
        newD=0
        TotalPopulation=countPopulation(ListePop)
        for j in range(0,len(ListePop)):
            nombreRepresentants=ListePop[j][0]
            if (nombreRepresentants>0):
                newCells=np.random.poisson(growthRate*nombreRepresentants*DT)
                newM=np.random.binomial(newCells,Ngenes*pm)
                if (newM>newCells):
                    newM=newCells
                
                newD=np.random.poisson(deathRate*nombreRepresentants*TotalPopulation*DT/KC)

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
        newH=np.random.poisson(pf*TotalPopulation*DT)
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

                        
                    (Off1,Off2)=fusionVect(Genotype1,Genotype2)
                    exist=0
                    count=0
                    while (count<len(ListePop) and exist==0):
                        if (convertBinaire(ListePop[count][1])-convertBinaire(Off1)==0):
                            ListePop[count][0]=ListePop[count][0]+1
                            exist=1
                        count=count+1
                    if (exist==0):
                        genotypesCounts=genotypesCounts+1


                        NewElement=[1,Off1,0,l*DT,ListePop[j][6],0,genotypesCounts]
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
                        NewElement=[1,Off2,0,l*DT,ListePop[j][6],0,genotypesCounts]
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

        
    f1First.write("%i\n" % (TotalCells[0]))
    f2First.write("%i\n" % (Number[0]))
    f3First.write("%i\n" % (Shanon[0]))
    f4First.write("%i\n" % (Simpson[0]))
    f5First.write("%i\n" % (Score[0]))

    f1First.close()
    f2First.close()
    f3First.close()
    f4First.close()
    f5First.close()

    f1.write("%i\n" % (TotalCells[NgenerationsMax-1]))
    f2.write("%i\n" % (Number[NgenerationsMax-1]))
    f3.write("%i\n" % (Shanon[NgenerationsMax-1]))
    f4.write("%i\n" % (Simpson[NgenerationsMax-1]))
    f5.write("%i\n" % (Score[NgenerationsMax-1]))

    #Closing file    
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    


s=sys.argv[1]


ModelRun(NgenerationsMax,DT,s,KC,pm,pf,diversity)

end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))