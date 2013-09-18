#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "fonks.cpp"

#define DEBUG 0
using namespace std;

/* Changelog
    0-0-2   2011/03/25 - correction erreurs (
                         psRank=sizeof(int) instead of size(char)
                         malloc cTransfert and psRank (+1);L128;
                         psL<valueNb+1 instead of psL<ValueNb;
                         process aValues[psL]==0: add valueNb++;psRank[valueNb]=psL;
    0-0-2   2010/06/14 - process network data written in a single file
    0-0-1   treat case c>0 and a=0
    0-0-0   first version : export to MCL and PAJEK NET format
*/

    int modelRank,netRank,psNb,cLimit;
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
    cTransfert=(char*) malloc(sizeof(char)*psNb);
    psRank=(int*) malloc(sizeof(int)*psNb);

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
    //process network file

    //open anti and corr files
    sprintf(inputFile,"c_m%03d_n%05d.4mat",modelRank,netRank);
    cReadFile=NULL;
    cReadFile=fopen(inputFile,"rb");
    if (cReadFile==NULL){
        printf("corr input file %s not opened\n",inputFile);
        exit(1);
    }
    sprintf(inputFile,"a_m%03d_n%05d.4mat",modelRank,netRank);
    aReadFile=NULL;
    aReadFile=fopen(inputFile,"rb");
    if (aReadFile==NULL){
        printf("anti input file %s not opened\n",inputFile);
        exit(1);
    }
    //Read C and A line by line

    for (lineL=0;lineL<psNb;lineL++){
        fread(aValues,sizeof(char),psNb,aReadFile);
        fread(cValues,sizeof(char),psNb,cReadFile);
        valueNb=-1;
        for (psL=0;psL<lineL;psL++){
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
                    valueNb++;
                    psRank[valueNb]=psL;
                    cTransfert[valueNb]=cValues[psL];
                }
            }
        }
        //process each value after x vs x
        for (psL=lineL+1;psL<psNb;psL++){
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
                    valueNb++;
                    psRank[valueNb]=psL;
                    cTransfert[valueNb]=cValues[psL];
                }
            }
        }
        //fwrite(transfert,sizeof(char),psNb,writeFile);
        for (psL=0;psL<valueNb+1;psL++){
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
    printf("end");
    fclose(writeFile);
    return 0;
}

int main(int argc, char *argv[])
{
    if (!DEBUG){
        int argNb;
        /*for(argNb=0;argNb<argc;argNb++){
            printf("argument %i: %s\n",argNb+1,argv[argNb]);
        }*/
        if(argc != 8){
            printf("needs modelRank netRank psNb  exportType cLimit cRatio  diffFlag\n");
            printf("/work/cinbell/exportnet 8 5 39485 0 1 500 MCL 20 1.0 1");
            return 0;
        }
        modelRank=atoi(argv[1]);
        netRank=atoi(argv[2]);
        psNb=atoi(argv[3]);
        copie_string(argv[4],exportType);
        cLimit=atoi(argv[5]);
        cRatio=atof(argv[6]);
        diffFlag=atoi(argv[7]);
        //printf("Doing exportation (modelRank: %d, netRank: %d\n",modelRank,netRank);
    }
    else{
        modelRank=11;
        netRank=42;
        psNb=22810;
        cLimit=30;
        cRatio=1.0;
        diffFlag=1;
        exportType[0]='M';
        exportType[1]='C';
        exportType[2]='L';
        exportType[3]='\0';
//      exportType[0]='P';
//      exportType[1]='A';
//      exportType[2]='J';
//      exportType[3]='\0';

        ;
    }
    printf("exportType: %s\n",exportType);
    if (strcmp(exportType,"MCL")==0|strcmp(exportType,"PAJ")==0){
        //printf("run export A\n");
        if (exportA()){
            //printf("error in exportA\n");
        }
    }
        //system("PAUSE");
    return 0;
}
