#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fonks.cpp"

#define DEBUG 0
using namespace std;

/* Changelog
    0-0-7   2012/12/14   add freqFlag (selection on % of changed in used comp)
    0-0-6   2012/12/10   add superior limit - correction rawFlag was not read
    0-0-5   2012/09/12   correct bug in TSN output
    0-0-4   2012/08/16   add local and distal server usage
    0-0-3   2012/02/21 - export subsets of probe sets
                         simplification of loops (avec if (psL!=lineL))
                         replace PAJEK by TSN (TensorNet)
    0-0-2   2011/03/25 - correction erreurs (
                         psRank=sizeof(int) istead of size(char)
                         malloc cTransfert and psRank (+1);L128;
                         psL<valueNb+1 inestead of psL<ValueNb;
                         process aValues[psL]==0: add valueNb++;psRank[valueNb]=psL;
    0-0-2   2010/06/14 - process network data written in a single file
    0-0-1   treat case c>0 and a=0
    0-0-0   first version : export to MCL and PAJEK NET format
*/

    int rankPos;
    unsigned int psNb,currPsNb,currPsRank,lineL,psL,distalFlag,listRank,Start;
    unsigned short modelRank,netRank,listL;
    char cInfLimit,cSupLimit,testValue,diffFlag,rawFlag,freqFlag;
    float cRatio;
    char exportType[256],inputFile[256],outputFile[256],fileName[256];
    char *cValues;
    char *aValues;
    char *fValues;
    char *cTransfert;
    unsigned int *psRanks;
    unsigned int *psRank;
    unsigned int lSize;
    unsigned int result;
    FILE *writeFile=NULL;
    FILE *writeFileRaw=NULL;
    FILE *aReadFile=NULL;
    FILE *cReadFile=NULL;
    FILE *fReadFile=NULL;

int exportA(){



    FILE *listFile;
    if (distalFlag==1){
        sprintf(fileName,"m%u_pslist%u.u32",modelRank,listRank);
    }
    else{
        sprintf(fileName,"/home/mbellis/array1/sosma/net/m%u/list/m%u_pslist%u.u32",modelRank,modelRank,listRank);
        printf("%s\n",fileName);
    }
    listFile = fopen ( fileName , "rb" );
    if (listFile==NULL) {fputs ("List file error",stderr); exit (1);}
    // obtain file size:
    fseek (listFile , 0 , SEEK_END);
    lSize = ftell (listFile);
    rewind (listFile);

    // allocate memory to contain the whole file:
    psRanks = (unsigned int*) malloc (sizeof(char)*lSize);
    if (psRanks == NULL) {fputs ("Memory error",stderr); exit (2);}

    // copy the file into the buffer:
    result = fread (psRanks,1,lSize,listFile);
    if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

    /* the whole file is now loaded in the memory buffer. */

    // terminate
    fclose (listFile);
    currPsNb=lSize/4;

    aValues=(char*) malloc(sizeof(char)*psNb);
    cValues=(char*) malloc(sizeof(char)*psNb);
    if (rawFlag or freqFlag){
        fValues=(char*) malloc(sizeof(char)*psNb);
    }
    //the x vs x value is not used
    cTransfert=(char*) malloc(sizeof(char)*currPsNb);
    psRank=(unsigned int*) malloc(sizeof(unsigned int)*currPsNb);


    //TotalNb=0;
    //open the output mcl file

    if (rawFlag){
        if (diffFlag){
            if (freqFlag){
                sprintf(outputFile,"m%d_n%d_c%dto%d_%s_%d_diff_raw_freq.txt",modelRank,netRank,cInfLimit,cSupLimit,exportType,listRank);
            }
            else{
                sprintf(outputFile,"m%d_n%d_c%dto%d_%s_%d_diff_raw.txt",modelRank,netRank,cInfLimit,cSupLimit,exportType,listRank);
            }
        }
        else{
            if (freqFlag){
                sprintf(outputFile,"m%d_n%d_c%dto%d_%s_%d_raw_freq.txt",modelRank,netRank,cInfLimit,cSupLimit,exportType,listRank);
            }
            else{
                sprintf(outputFile,"m%d_n%d_c%dto%d_%s_%d_raw.txt",modelRank,netRank,cInfLimit,cSupLimit,exportType,listRank);
            }
        }
    }
    else{
        if (diffFlag){
            if (freqFlag){
                sprintf(outputFile,"m%d_n%d_c%dto%d_%s_%d_diff_freq.txt",modelRank,netRank,cInfLimit,cSupLimit,exportType,listRank);
            }
            else{
                sprintf(outputFile,"m%d_n%d_c%dto%d_%s_%d_diff.txt",modelRank,netRank,cInfLimit,cSupLimit,exportType,listRank);
            }
        }
        else{
            if (freqFlag){
                sprintf(outputFile,"m%d_n%d_c%dto%d_%s_%d_freq.txt",modelRank,netRank,cInfLimit,cSupLimit,exportType,listRank);
            }
            else{
                sprintf(outputFile,"m%d_n%d_c%dto%d_%s_%d.txt",modelRank,netRank,cInfLimit,cSupLimit,exportType,listRank);
            }
        }
    }

    writeFile=NULL;
    writeFile=fopen(outputFile,"w");
    if (writeFile==NULL){
        printf("output file %s not opened",outputFile);
        exit(1);
    }
    //process network file

    //open anti and corr files
    if (distalFlag==1){
        sprintf(inputFile,"c_m%d_n%d.4mat",modelRank,netRank);
    }
    else{
        sprintf(inputFile,"/home/mbellis/array1/sosma/net/m%u/n%u/c_m%d_n%d.4mat",modelRank,netRank,modelRank,netRank);
        printf("%s\n",inputFile);
    }
    cReadFile=NULL;
    cReadFile=fopen(inputFile,"rb");
    if (cReadFile==NULL){
        printf("corr input file %s not opened\n",inputFile);
        exit(1);
    }
    if (distalFlag==1){
        sprintf(inputFile,"a_m%d_n%d.4mat",modelRank,netRank);
        }
    else{
        sprintf(inputFile,"/home/mbellis/array1/sosma/net/m%u/n%u/a_m%d_n%d.4mat",modelRank,netRank,modelRank,netRank);
        printf("%s\n",inputFile);
    }
    aReadFile=NULL;
    aReadFile=fopen(inputFile,"rb");
    if (aReadFile==NULL){
        printf("anti input file %s not opened\n",inputFile);
        exit(1);
    }
    if (rawFlag or freqFlag){
        sprintf(inputFile,"f_m%d_n%d.4mat",modelRank,netRank);
        fReadFile=NULL;
        fReadFile=fopen(inputFile,"rb");
        if (fReadFile==NULL){
            printf("freq input file %s not opened\n",inputFile);
            exit(1);
        }
    }
    //Read C and A line by line
    for (lineL=0;lineL<currPsNb;lineL++){
        if (strcmp(exportType,"MCL")==0){
            Start=0;
        }
        else if (strcmp(exportType,"TSN")==0){
            Start=lineL+1;
        }

    //for (lineL=currPsNb-1;lineL<currPsNb;lineL++){
    //for (lineL=0;lineL<1;lineL++){
        currPsRank=psRanks[lineL]-1;
        fseek(aReadFile,sizeof(char)*currPsRank*psNb,SEEK_SET);
        fseek(cReadFile,sizeof(char)*currPsRank*psNb,SEEK_SET);
        fread(aValues,sizeof(char),psNb,aReadFile);
        fread(cValues,sizeof(char),psNb,cReadFile);
        if (rawFlag or freqFlag){
            fseek(fReadFile,sizeof(char)*currPsRank*psNb,SEEK_SET);
            fread(fValues,sizeof(char),psNb,fReadFile);
        }
        rankPos=-1;
        for (psL=Start;psL<currPsNb;psL++){
            currPsRank=psRanks[psL]-1;
            if (psL!=lineL){
                if (rawFlag){
                    if (diffFlag==0){
                        testValue=(cValues[currPsRank]*fValues[currPsRank])/100;
                    }
                    else {
                        testValue=((cValues[currPsRank]-aValues[currPsRank])*fValues[currPsRank])/100;
                    }
                }
                else{
                    if (freqFlag){
                        testValue=fValues[currPsRank];
                    }
                    else if (diffFlag==0){
                        testValue=cValues[currPsRank];
                    }
                    else {
                        testValue=cValues[currPsRank]-aValues[currPsRank];
                    }
                }
                if ((testValue>cInfLimit)&(testValue<=cSupLimit)){
                    if (aValues[currPsRank]>0){
                        if ((float)cValues[currPsRank]/(float)aValues[currPsRank]>=cRatio){
                            rankPos++;
                            if (strcmp(exportType,"MCL")==0){
                                psRank[rankPos]=psRanks[psL];
                            }
                            else if (strcmp(exportType,"TSN")==0){
                                psRank[rankPos]=psL;
                            }
                            cTransfert[rankPos]=testValue;
                        }
                    }
                    else{
                        rankPos++;
                        if (strcmp(exportType,"MCL")==0){
                            psRank[rankPos]=psRanks[psL];
                        }
                        else if (strcmp(exportType,"TSN")==0){
                            psRank[rankPos]=psL;
                        }
                        cTransfert[rankPos]=testValue;
                    }
                }
            }
        }
        //fwrite(transfert,sizeof(char),currPsNb,writeFile);
        for (psL=0;psL<(unsigned int)(rankPos+1);psL++){
            if (strcmp(exportType,"MCL")==0){
                fprintf(writeFile,"%d %d %d\n",psRanks[lineL],psRank[psL],cTransfert[psL]);
            }
            else if (strcmp(exportType,"TSN")==0){
                fprintf(writeFile,"%d\t%d\t%0.2f\n",lineL+1,psRank[psL]+1,(float)cTransfert[psL]/100.0);
            }
        }
    }
    fclose(writeFile);
    fclose(aReadFile);
    fclose(cReadFile);
    free(aValues);
    free(cValues);
    if (rawFlag){
        fclose(fReadFile);
        free(fValues);
    }
    free(cTransfert);
    free (psRanks);

    return 0;
}

int main(int argc, char *argv[])
{
    if (!DEBUG){
        /*int argNb;
        for(argNb=0;argNb<argc;argNb++){
            printf("argument %i: %s\n",argNb+1,argv[argNb]);
        }*/
        if(argc != 13){
            printf("needs modelRank netRank psNb  exportType cInfLimit cSupLimit cRatio  diffFlag listRank rawFlag distalFlag freqFlag\n");
            printf("/work/cinbell/exportnet 27 22 22690 MCL 30 70 1.0 1 5 0 1\n");
            printf("/work/cinbell/exportnet 8 24 45101 MCL 30 70 1.0 1 5 1 0\n");
            return 0;
        }
        modelRank=atoi(argv[1]);
        netRank=atoi(argv[2]);
        psNb=atol(argv[3]);
        copie_string(argv[4],exportType);
        cInfLimit=atoi(argv[5]);
        cSupLimit=atoi(argv[6]);
        cRatio=atof(argv[7]);
        diffFlag=atoi(argv[8]);
        listRank=atoi(argv[9]);
        rawFlag=atoi(argv[10]);
        distalFlag=atoi(argv[11]);
        freqFlag=atoi(argv[11]);

    }
    else{
        modelRank=27;
        netRank=22;
        psNb=22690;
        cInfLimit=30;
        cRatio=1.0;
        diffFlag=1;
        exportType[0]='M';
        exportType[1]='C';
        exportType[2]='L';
        exportType[3]='\0';
        listRank=4;
        distalFlag=1;
        //exportnet.exe 27 22 22690 30 1 1 MCL
//      exportType[0]='P';
//      exportType[1]='A';
//      exportType[2]='J';
//      exportType[3]='\0';

        ;
    }
    if ((strcmp(exportType,"MCL")==0)|(strcmp(exportType,"TSN")==0)){
        //printf("run export A\n");
        if (exportA()){
            //printf("error in exportA\n");
        }
    }
        //system("PAUSE");
    return 0;
}
