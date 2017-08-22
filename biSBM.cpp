// biSBM v1.2
//
//Daniel Larremore
//Harvard School of Public Health
//July 29, 2014
//http://danlarremore.com/bipartiteSBM
//daniel.larremore@gmail.com
//
//biSBM - a method for community detection in bipartite networks, based on the publication:
//"Efficiently inferring community structure in bipartite networks"
//Daniel B. Larremore, Aaron Clauset, and Abigail Z. Jacobs
//Physical Review E 90(1), 012805 (2014).
//http://danlarremore.com/pdf/2014_LCJ_EfficientlyInferringCommunityStructureInBipartiteNetworks_PRE.pdf
//
// Please do not distribute without contacting the author above at daniel.larremore@gmail.com
// If a bug is located within the code, please contact the author, to correct the official version!
//
// This code is based on code written by Brian Karrer for the stochastic block model, http://arxiv.org/abs/1104.3590
// You can download that code at http://www-personal.umich.edu/~mejn/dcbm/
//

// ***** TO COMPILE *****
// g++ -O3 -Wall -g -pedantic -o biSBM biSBM.cpp
// **********************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>
#include <vector>
#include "biSBM.h"

using namespace std;
string trueCommsName;

int main(int argc, char *argv[])
{
    srandom(time(NULL));
    // Demand correct usage
    //printf("%i\n",argc);
    if (argc==1)
    {
        //                 1           2               3            4    5          5             6         7                      8
        printf("Usage: %s <edgeList> <vertexTypes> <outputFOLDER> <Ka>  <Kb> <isDegreeCorrect> <KL steps> <optionalTrueComms> <initialize@truecomms?>\n", argv[0]);
        return 1;
    }
    // Read in argv from function call
    string edgeListName = argv[1];
    string vertexTypesName = argv[2];
    string folderName = argv[3];
    isDegreeCorrect = 0;
    
    // Create the Comms vector of vectors
    // Comms[i] is a vector of community numbers corresponding to type i
    int counter = 0;
        vector<int> commlist;
        for (unsigned int q=0; q<strtol(argv[4],NULL,10); ++q) {
            commlist.push_back(counter);
            counter++;
        }
        Comms.push_back(commlist);
    commlist.clear();
    for (unsigned int q=0; q<strtol(argv[5],NULL,10); ++q) {
        commlist.push_back(counter);
        counter++;
    }
    Comms.push_back(commlist);
    
    MaxComms = counter;
    cout << "*****" << " Communities by type/part" << " *****" << endl;
    for (unsigned int i=0; i<Comms.size(); ++i) {
        for (unsigned int j=0; j<Comms[i].size(); ++j) {
            printf("%i,",Comms[i][j]);
        }
        printf("\n");
    }

    
    if (argc>=7)
    {
        isDegreeCorrect = strtol(argv[6],NULL,10);
        if (isDegreeCorrect==1)
        {
            printf("Using biSBM with degree correction.\n");
        }
        else if (isDegreeCorrect==0)
        {
            printf("Using biSBM without degree correction.\n");
        }
        else
        {
            printf("Unknown parameter for isDegreeCorrect. Please use 1 or 0.\n");
            return 1;
        }
        
        if (argc>=8)
        {
            KLPerNetwork = strtol(argv[7],NULL,10);
        }
        
        if (argc>=9)
        {
            trueCommsName = argv[8];
            TrueCommsAvailable = 1;
            printf("TRUE COMMS FROM: %s\n",trueCommsName.c_str());
            if (argc>=9)
            {
                InitializationOption = strtol(argv[9],NULL,10);
                printf("INITIALIZING AT TRUE COMMS.\n");
            }
            else
            {
                printf("INITIALIZING AT RANDOM COMMS: %s\n",trueCommsName.c_str());
            }
        }
    }
    else
    {
        printf("Using SBM.\n");
    }
    
    printf("edges: \t%s\n",edgeListName.c_str());
    printf("types: \t%s\n",vertexTypesName.c_str());
    printf("K: \t%s,%s\n",argv[4],argv[5]);
    printf("KL steps\t%d\n",KLPerNetwork);
    printf("output: \t%s\n",folderName.c_str());
    mkdir(folderName.c_str(),0777);
    
    unsigned int i, j, k;
    
    // Read in the network
    GetTheVertexTypes(vertexTypesName);
    GetTheNetworkEdges(edgeListName);
    
    HighestScore = -numeric_limits<double>::max( );
    VIValue = 0;
    NMIValue = 0;
    time_t startTime = time(NULL);
    ofstream logfile;
    logfile.open((folderName + "/logfile").c_str());
    for(j=0; j < KLPerNetwork; j++)
    {
        RunKL();
        //KL,dt,L:
        logfile << j+1 << "," << difftime(time(NULL),startTime) << "," << max(MaxScore,HighestScore) << endl;
        printf(">%i,%f,%f\n",j+1,difftime(time(NULL),startTime), max(MaxScore,HighestScore));
        
        if(MaxScore >= HighestScore)
        {
            HighestScore = MaxScore;
            if(TrueCommsAvailable == 1)
            {
                VIValue = ComputeVI();
                NMIValue = ComputeNMI();
                cout << "VI Value: " << VIValue << " NMI Value: " << NMIValue << endl;
            }
            for(i=0; i < MaxComms; i++)
            {
                SavedCommVertices[i] = BestCommVertices[i];
                SavedCommStubs[i] = BestCommStubs[i];
                for(k=0; k < MaxComms; k++)
                    SavedEdgeMatrix[i][k] = BestEdgeMatrix[i][k];
            }
            for(i=0; i < Nodes; i++)
                SavedState[i] = BestState[i];
        }
    }
    logfile.close();
    
    // because PrintResults are written for best values we copy them
    // back over from the saved values which are the best ones.
    for(i=0; i < MaxComms; i++)
    {
        BestCommVertices[i] = SavedCommVertices[i];
        BestCommStubs[i] = SavedCommStubs[i];
        for(k=0; k < MaxComms; k++)
            BestEdgeMatrix[i][k] = SavedEdgeMatrix[i][k];
    }
    for(i=0; i < Nodes; i++)
        BestState[i] = SavedState[i];
    cout << "Final Score: " << ComputeInitialScore() << endl;
    
    PrintResults(folderName);
    
    
    freegraph();
    
    return 0;
}
