#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "fonks.cpp"

#define DEBUG 0
using namespace std;

/* Changelog
    0-0-1   treat case c>0 and a=0
    0-0-0   first version : export to MCL and PAJEK NET format
*/

    int modelRank,netRank,psNb,cLimit,firstBloc,lastBloc,blocSize,lastBlocSize;
    int diffFlag;
    float cRatio;
    int blocL,lineL,psL,offSet,currBlocSize,valueNb;
    char exportType[256];
    char inputFile[256],outputFile[256];
    char *cValues;
    char *aValues;
    char *cTransfert;
    int *psRank;
    FILE *writeFile=NULL;
    FILE *aReadFile=NULL;
    FILE *cReadFile=NULL;

int exportA(){

    aValues=(char*) malloc(sizeof(char)*psNb);
    cValues=(char*) malloc(sizeof(char)*psNb);
    //the x vs x value is not used
    cTransfert=(char*) malloc(sizeof(char)*(psNb-1));
    psRank=(int*) malloc(sizeof(char)*(psNb-1));

    //TotalNb=0;
    //open the output mcl file
    if (strcmp(exportType,"MCL")==0){
        //printf("starting MCL\n");
        sprintf(outputFile,"m%03d_n%05d_c%d_mcl.txt",modelRank,netRank,cLimit);
    }
    else if (strcmp(exportType,"PAJ")==0){
        //printf("starting PAJ\n");
        sprintf(outputFile,"m%03d_n%05d_c%d_paj.txt",modelRank,netRank,cLimit);
    }
    writeFile=NULL;
    writeFile=fopen(outputFile,"w");
    if (writeFile==NULL){
        printf("output file %s not opened",outputFile);
        exit(1);
    }
    if (strcmp(exportType,"PAJ")==0){
        fprintf(writeFile,"*vertices %d\n",psNb);
        for (psL=0;psL<psNb;psL++){
            fprintf(writeFile,"%d \"%d\"\n",psL+1,psL+1);
        }
        fprintf(writeFile,"*edges\n");
    }
    //process each bloc
    for (blocL=firstBloc;blocL<lastBloc;blocL++){
        offSet=blocSize*blocL;
        //open anti and corr files
        sprintf(inputFile,"cm_m%03d_n%05d_%03d.4mat",modelRank,netRank,blocL);
        cReadFile=NULL;
        cReadFile=fopen(inputFile,"rb");
        if (cReadFile==NULL){
            printf("corr input file %s not opened",inputFile);
            exit(1);
        }
        sprintf(inputFile,"am_m%03d_n%05d_%03d.4mat",modelRank,netRank,blocL);
        aReadFile=NULL;
        aReadFile=fopen(inputFile,"rb");
        if (aReadFile==NULL){
            printf("anti input file %s not opened",inputFile);
            exit(1);
        }
        //Read C and A line by line
        if (blocL<lastBloc-1){
            currBlocSize=blocSize;
        }
        else{
            currBlocSize=lastBlocSize;
        }
        for (lineL=0;lineL<currBlocSize;lineL++){
            fread(aValues,sizeof(char),psNb,aReadFile);
            fread(cValues,sizeof(char),psNb,cReadFile);
//            for (psL=0;psL<50;psL++){
//                printf("%d %d %d\n",psL,aValues[psL],cValues[psL]);
//            }
//            exit(EXIT_FAILURE);
            //process each value before x vs x
            valueNb=-1;
            for (psL=0;psL<blocSize*blocL+lineL;psL++){
                if (cValues[psL]>cLimit){
                    if (aValues[psL]>0){
                        if ((float)cValues[psL]/aValues[psL]>=cRatio){
                            valueNb++;
                            psRank[valueNb]=psL;
                            if (diffFlag==1){
                                cTransfert[valueNb]=cValues[psL]-aValues[psL];
                            }
                            else{
                                cTransfert[valueNb]=cValues[psL];
                            }
                        }
                    }
                    else{
                        cTransfert[valueNb]=cValues[psL];
                    }
                }
            }
            //process each value after x vs x
            for (psL=blocSize*blocL+lineL+1;psL<psNb;psL++){
                if (cValues[psL]>cLimit){
                    if (aValues[psL]>0){
                        if ((float)cValues[psL]/aValues[psL]>=cRatio){
                            valueNb++;
                            psRank[valueNb]=psL;
                            if (diffFlag==1){
                                cTransfert[valueNb]=cValues[psL]-aValues[psL];
                            }
                            else{
                                cTransfert[valueNb]=cValues[psL];
                            }
                        }
                    }
                    else{
                        cTransfert[valueNb]=cValues[psL];
                    }
                }
            }
            //fwrite(transfert,sizeof(char),psNb,writeFile);
            for (psL=0;psL<valueNb;psL++){
                if (strcmp(exportType,"MCL")==0){
                    fprintf(writeFile,"%d %d %d\n",lineL+offSet+1,psRank[psL]+1,cTransfert[psL]);
                }
                else if (strcmp(exportType,"PAJ")==0){
                    fprintf(writeFile,"%d %d %0.2f\n",lineL+offSet+1,psRank[psL]+1,cTransfert[psL]/100.0);
                }
            }
        }
        fclose(aReadFile);
        fclose(cReadFile);
    }
    fclose(writeFile);
    return 0;
}

int main(int argc, char *argv[])
{
    if (!DEBUG){
        int argNb;
        for(argNb=0;argNb<argc;argNb++){
            printf("argument %i: %s\n",argNb+1,argv[argNb]);
        }
        if(argc != 11){
            printf("needs modelRank netRank psNb firstBloc lastBloc blocSize exportType cLimit cRatio  diffFlag\n");
            printf("/work/cinbell/exportnet 8 5 39485 0 1 500 MCL 20 1.0 1");
            return 0;
        }
        modelRank=atoi(argv[1]);
        netRank=atoi(argv[2]);
        psNb=atoi(argv[3]);
        firstBloc=atoi(argv[4]);
        lastBloc=atoi(argv[5]);
        blocSize=atoi(argv[6]);
        copie_string(argv[7],exportType);
        cLimit=atoi(argv[8]);
        cRatio=atof(argv[9]);
        diffFlag=atoi(argv[10]);
        printf("Doing exportation (modelRank: %d, netRank: %d, firstBloc: %d and lastBloc: %d, bloc size: %d \n",modelRank,netRank,firstBloc,lastBloc,blocSize);
    }
    else{
        modelRank=8;
        netRank=5;
        psNb=39485;
        cLimit=20;
        cRatio=1.0;
        firstBloc=0;
        lastBloc=1;
        blocSize=500;
        diffFlag=1;
//        exportType[0]='M';
//        exportType[1]='C';
//        exportType[2]='L';
//        exportType[3]='\0';
        exportType[0]='P';
        exportType[1]='A';
        exportType[2]='J';
        exportType[3]='\0';

        ;
    }
    lastBlocSize=psNb%blocSize;
    printf("exportType: %s - lastBlocSize: %d\n",exportType,lastBlocSize);
    if (strcmp(exportType,"MCL")==0|strcmp(exportType,"PAJ")==0){
        printf("run export A\n");
        if (exportA()){
            printf("error in exportA\n");
        }
    }
        //system("PAUSE");
    return 0;
}
