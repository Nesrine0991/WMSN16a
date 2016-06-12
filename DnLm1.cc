#include "DnLm1.h"
#include <cassert>
#include <cstdint>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include "NetwControlInfo.h"
#include "DnLm1Pkt_m.h"
#include "SimTracer.h"
#include "DownHillSimplex.h"
#include <string>
#include <fstream>
#include <chrono>
#include <array>
#include <vector>
#include <iostream>

Define_Module(DnLm1);

void DnLm1::initialize(int stage)
{
    BaseNetwLayer::initialize(stage);

    if(stage == 1) {
        sinkAddress =9;
        stats = par("stats");
        trace = par("trace");
        debug = par("debug");
        floodSeqNumber = 0;
        nbDataPacketsForwarded = 0;
        nbDataPacketsReceived = 0;
        nbDataPacketsSent = 0;
        nbHops = 0;
        i=0;
        l=0;
        N=10;
        route=0;
        link=22;
        V=9;
         h=1;
         k=0;
         alpha=0.5;
         beta=1.3*pow(10.0,-8.0);
         sigma2=3500;
         gam=55.54;
         np=4;
         cr=0.5;
         bi=5.0;
         w=0.15;
         repV=0;
         Dh=100;
         deltap=0.2;
         deltar=0.2;
         deltax=0.2;
         epsilon=1E-10;
         count=1;
         verification=0;
         sumqiii=0;
         somun=0;
         sumxhl=0;
         int ailorg[10][22]={
              1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              -1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, -1, 0, -1, 0, -1, 0, 0, -1, 0, 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 1, 1, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, -1, 0, 1,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, -1, -1
        };

        ail = (int**)malloc (10*sizeof (int*));
        for (i=0;i<10;i++){
            ail [i]= (int*)malloc ((22)*sizeof (int));
            for(j=0;j<22;j++){
                ail[i][j]=ailorg[i][j];
            }
        }
        double wlorg[22]={0.015982,0.016098,0.018345,0.027314,0.026921,0.029343,0.015215,0.016234,0.019233,0.017234,0.0128671, 0.0207475, 0.0288934, 0.0251723, 0.020437, 0.0121182,0.0223527, 0.0129935, 0.0222419, 0.0144464, 0.0223387, 0.0177298};
        wl = (double*)malloc (22*sizeof (double));
        for (l=0;l<link;l++){
            if(ail[myNetwAddr][l]==0)
                wl[l]=0;
            else
                wl[l]=wlorg[l];
         }
        double qiorg[10]={0.19982,0.16098,0.18345,0.15314,0.17923,0.17343,0.19215,0.18834,0.18233,0.19234};
        double psorg[10]={3.05,3.05,3.05,3.05,3.05,3.05,3.05,3.05,3.05,3.05};
        double lamdaorg[10]={0.209763,0.218569,0.243038,0.268853,0.220553,0.271589,0.208977,0.26945,0.184731,0.224713};
        double vhorg[10]={0.209763,0.218569,0.243038,0.268853,0.220553,0.271589,0.208977,0.26945,0.184731,0.204713};
        double rorg[10]={0.0304375,0.0304375,0.0304375,0.0304375,0.0304375,0.0304375,0.0304375,0.0304375,0.0304375,0.0304375};
        double degretab[9]={3,3,3,4,2,3,1,2,1};

        qi = (double*)malloc (10*sizeof (double));
        ps = (double*)malloc (10*sizeof (double));
        lamda = (double*)malloc (10*sizeof (double));
        vh = (double*)malloc (10*sizeof (double));
        r = (double*)malloc (10*sizeof (double));
        for (i=0;i<N;i++){
                qi[i]=qiorg[i];
                ps[i]=psorg[i];
                lamda[i]=lamdaorg[i];
                vh[i]=vhorg[i];
                r[i]=rorg[i];
                RhVoisin[i]=0;
            }
        l=0;


        double uhiorg[9][10]={
           0.209763, 0.218569, 0.243038, 0.268853, 0.220553, 0.271589, 0.208977, 0.26945, 0.184731, 0.224713,
           0.229179, 0.176876, 0.187517, 0.159507, 0.278355, 0.111343, 0.292733, 0.154531, 0.176688, 0.195533,
           0.258345, 0.262434, 0.205779, 0.195995, 0.213609, 0.178557, 0.285119, 0.267216, 0.114207, 0.167479,
           0.117426, 0.229634, 0.104044, 0.173648, 0.266524, 0.291431, 0.255631, 0.12807, 0.274002, 0.274017,
           0.295724, 0.194722, 0.259832, 0.260182, 0.192296, 0.204095, 0.256106, 0.235776, 0.123655, 0.244127,
           0.227984, 0.216404, 0.128671, 0.207475, 0.288934, 0.251723, 0.20437, 0.121182, 0.182932, 0.19472,
           0.152911, 0.137266, 0.254847, 0.247384, 0.19123, 0.14331, 0.213687, 0.127044, 0.103758, 0.164828,
           0.223527, 0.129935, 0.222419, 0.144464, 0.223387, 0.177298, 0.28875, 0.28052, 0.236364, 0.18999,
           0.171902, 0.222613, 0.187406, 0.28047, 0.239526, 0.119856, 0.112045, 0.293962, 0.233353, 0.230628,
         };

        uhi = (double**)malloc (9*sizeof (double*));
        for (i=0;i<9;i++){
            uhi [i]= (double*)malloc ((10)*sizeof (double));
            for(j=0;j<10;j++){
                uhi[i][j]= uhiorg[i][j];
            }
        }


        for(int i=0;i<N;i++){
          for(int j=0;j<link;j++){
              if(ail[i][j]==1) {
                 ailp[i][j]=1;
                 ailm[i][j]=0;
              }
              else
                 if(ail[i][j]==-1){
                    ailm[i][j]=1;
                    ailp[i][j]=0;
                 }
                 else{
                    ailm[i][j]=0;
                    ailp[i][j]=0;
                 }
           }
         }
        for(i=0;i<16;i++){
            tabVoisin[i]=-1;
            tablien[i]=-1;
        }
        for(i=0;i<9;i++){
            for(j=0;j<22;j++){
                xhlVoisin[i][j]=0.0;
            }

        }
        for(i=0;i<N;i++){
            qiV[i]=0;
            qiVV[i]=0;
        }
        for(i=0;i<link;i++)
            xhl[i]=0;
        

        double dist[22]={20.1, 14.14, 26.92,
                14.21, 15.52, 18.35,
                25.71, 13.04, 22.82,
                13.6, 25.08, 21.59, 22.47,
                18.38, 27.5,
                13.41, 28.6, 19.03,
                10.39,
                17.52, 16.64,
                13.6};

        for (i=0;i<link;i++){
            cls[i]=alpha+beta*pow(dist[i],np);
         }


        if(myNetwAddr==0){ tabVoisin[0]=2; tabVoisin[1]=20.1;tabVoisin[2]=3; tabVoisin[3]=14.14;tabVoisin[4]=4; tabVoisin[5]=26.92;}
        if(myNetwAddr==1){ tabVoisin[0]=3; tabVoisin[1]=14.21; tabVoisin[2]=4; tabVoisin[3]=15.52; tabVoisin[4]=5; tabVoisin[5]=18.35;}
        if(myNetwAddr==2){ tabVoisin[0]=4; tabVoisin[1]=25.71; tabVoisin[2]=6; tabVoisin[3]=13.04; tabVoisin[4]=7; tabVoisin[5]=22.82; tabVoisin[6]=0; tabVoisin[7]=22.36;}
        if(myNetwAddr==3){ tabVoisin[0]=4; tabVoisin[1]=13.6; tabVoisin[2]=5; tabVoisin[3]=25.08; tabVoisin[4]=6; tabVoisin[5]=21.59;tabVoisin[6]=7; tabVoisin[7]=22.47 ;tabVoisin[8]=0; tabVoisin[9]=14.14; tabVoisin[10]=1; tabVoisin[11]=16.28;}
        if(myNetwAddr==4){ tabVoisin[0]=8; tabVoisin[1]=18.38; tabVoisin[2]=9; tabVoisin[3]=27.5;tabVoisin[4]=1; tabVoisin[5]=15.52;tabVoisin[6]=2; tabVoisin[7]=25.7;tabVoisin[8]=3; tabVoisin[9]=13.6;tabVoisin[10]=5; tabVoisin[11]=13.41;tabVoisin[12]=0; tabVoisin[13]=26.92;}
        if(myNetwAddr==5){ tabVoisin[0]=4; tabVoisin[1]=13.41; tabVoisin[2]=7; tabVoisin[3]=28.6; tabVoisin[4]=8; tabVoisin[5]=19.03;tabVoisin[6]=1; tabVoisin[7]=18.35;tabVoisin[8]=3; tabVoisin[9]=25.08;}
        if(myNetwAddr==6){ tabVoisin[0]=9; tabVoisin[1]=10.39;tabVoisin[2]=2; tabVoisin[3]=13.04;tabVoisin[4]=3; tabVoisin[7]=21.59;}
        if(myNetwAddr==7){ tabVoisin[0]=8; tabVoisin[1]=17.52; tabVoisin[2]=9; tabVoisin[3]=16.64;tabVoisin[4]=2; tabVoisin[5]=22.82;tabVoisin[6]=3; tabVoisin[7]=22.47;tabVoisin[8]=5; tabVoisin[9]=28.6;}
        if(myNetwAddr==8){ tabVoisin[0]=9; tabVoisin[1]=13.6;tabVoisin[2]=4; tabVoisin[3]=18.38;tabVoisin[4]=5; tabVoisin[5]=9.03;tabVoisin[6]=7; tabVoisin[7]=17.52;}
        if(myNetwAddr==9){ tabVoisin[0]=4; tabVoisin[1]=27.5;tabVoisin[2]=6; tabVoisin[3]=10.39;tabVoisin[4]=7; tabVoisin[5]=16.64;tabVoisin[6]=8; tabVoisin[7]=13.6;}


        if(myNetwAddr==0){ xhl[0]=0.34414; xhl[1]=0.37050; xhl[2]=0.34414;}
        if(myNetwAddr==1){ xhl[3]=0.34414; xhl[4]=0.37050; xhl[5]=0.34414;}
        if(myNetwAddr==2){ xhl[6]=0.35195; xhl[7]=0.34414; xhl[8]=0.36367; xhl[0]=0.34414;}
        if(myNetwAddr==3){ xhl[9]=0.34414; xhl[10]=0.35195; xhl[11]=0.35155; xhl[12]=0.36367; xhl[1]=0.37050;}
        if(myNetwAddr==4){ xhl[13]=0.37050; xhl[14]=0.36367; xhl[2]=0.34414; xhl[4]=0.37050; xhl[6]=0.35195; xhl[9]=0.34414; xhl[15]=0.35195;}
        if(myNetwAddr==5){ xhl[15]=0.35195; xhl[16]=0.36367; xhl[17]=0.37050; xhl[5]=0.34414; xhl[10]=0.35195;}
        if(myNetwAddr==6){ xhl[18]=0.36367; xhl[7]=0.34414; xhl[11]=0.35155;}
        if(myNetwAddr==7){ xhl[19]=0.37050; xhl[20]=0.36367; xhl[8]=0.36367; xhl[12]=0.36367; xhl[16]=0.36367;}
        if(myNetwAddr==8){ xhl[21]=0.34414; xhl[13]=0.37050; xhl[17]=0.37050; xhl[19]=0.37050;}
        if(myNetwAddr==9){ xhl[14]=0.36367; xhl[18]=0.36367; xhl[20]=0.36367; xhl[21]=0.34414;}

        if(myNetwAddr==0){ xhlVoisin[2][0]=0.34414; xhlVoisin[3][1]=0.37050; xhlVoisin[4][2]=0.34414;}
        if(myNetwAddr==1){ xhlVoisin[3][3]=0.34414; xhlVoisin[4][4]=0.37050; xhlVoisin[5][5]=0.34414;}
        if(myNetwAddr==2){ xhlVoisin[4][6]=0.35195; xhlVoisin[6][7]=0.34414; xhlVoisin[7][8]=0.36367; xhlVoisin[0][0]=0.34414;}
        if(myNetwAddr==3){ xhlVoisin[4][9]=0.34414; xhlVoisin[5][10]=0.35195; xhlVoisin[6][11]=0.35155; xhlVoisin[7][12]=0.36367; xhlVoisin[0][1]=0.37050; xhlVoisin[1][3]=0.34414;;}
        if(myNetwAddr==4){ xhlVoisin[8][13]=0.37050; xhlVoisin[9][14]=0.36367; xhlVoisin[0][2]=0.34414; xhlVoisin[1][4]=0.37050; xhlVoisin[2][6]=0.35195; xhlVoisin[3][9]=0.34414; xhlVoisin[5][15]=0.35195;}
        if(myNetwAddr==5){ xhlVoisin[4][15]=0.35195; xhlVoisin[7][16]=0.36367; xhlVoisin[8][17]=0.37050; xhlVoisin[1][5]=0.34414; xhlVoisin[3][10]=0.35195;}
        if(myNetwAddr==6){ xhlVoisin[9][18]=0.36367; xhlVoisin[2][7]=0.34414; xhlVoisin[3][11]=0.35155;}
        if(myNetwAddr==7){ xhlVoisin[8][19]=0.37050; xhlVoisin[9][20]=0.36367; xhlVoisin[2][8]=0.36367; xhlVoisin[3][12]=0.36367; xhlVoisin[5][16]=0.36367;}
        if(myNetwAddr==8){ xhlVoisin[9][21]=0.34414; xhlVoisin[4][13]=0.37050; xhlVoisin[5][17]=0.37050; xhlVoisin[7][19]=0.37050;}
        if(myNetwAddr==9){ xhlVoisin[8][21]=0.34414; xhlVoisin[7][20]=0.36367; xhlVoisin[6][18]=0.36367; xhlVoisin[4][14]=0.36367;}

        if(myNetwAddr==0){ tablien[0]=2; tablien[1]=0;tablien[2]=3; tablien[3]=1;tablien[4]=4; tablien[5]=2;}
        if(myNetwAddr==1){ tablien[0]=3; tablien[1]=3; tablien[2]=4; tablien[3]=4; tablien[4]=5; tablien[5]=5;}
        if(myNetwAddr==2){ tablien[0]=4; tablien[1]=6; tablien[2]=6; tablien[3]=7; tablien[4]=7; tablien[5]=8; tablien[6]=0; tablien[7]=0;}
        if(myNetwAddr==3){ tablien[0]=4; tablien[1]=9; tablien[2]=5; tablien[3]=10; tablien[4]=6; tablien[5]=11;tablien[6]=7; tablien[7]=12 ;tablien[8]=0; tablien[9]=1; tablien[10]=1; tablien[11]=3;}
        if(myNetwAddr==4){ tablien[0]=8; tablien[1]=13; tablien[2]=9; tablien[3]=14;tablien[4]=1; tablien[5]=4;tablien[6]=2; tablien[7]=6;tablien[8]=3; tablien[9]=9;tablien[10]=5; tablien[11]=15;tablien[12]=0; tablien[13]=2;}
        if(myNetwAddr==5){ tablien[0]=4; tablien[1]=15; tablien[2]=7; tablien[3]=16; tablien[4]=8; tablien[5]=17;tablien[6]=1; tablien[7]=5;tablien[8]=3; tablien[9]=10;}
        if(myNetwAddr==6){ tablien[0]=9; tablien[1]=18;tablien[2]=2; tablien[3]=7;tablien[4]=3; tablien[7]=11;}
        if(myNetwAddr==7){ tablien[0]=8; tablien[1]=19; tablien[2]=9; tablien[3]=20;tablien[4]=2; tablien[5]=8;tablien[6]=3; tablien[7]=12;tablien[8]=5; tablien[9]=16;}
        if(myNetwAddr==8){ tablien[0]=9; tablien[1]=21;tablien[2]=4; tablien[3]=13;tablien[4]=5; tablien[5]=17;tablien[6]=7; tablien[7]=19;}
        if(myNetwAddr==9){ tablien[0]=4; tablien[1]=14;tablien[2]=6; tablien[3]=18;tablien[4]=7; tablien[5]=20;tablien[6]=8; tablien[7]=21;}

        for( h=0;h<N;h++){
          i=0; rep=0;
          while (tabVoisin[i]!=-1 ){
              if (tabVoisin[i]==h)
                  rep=1;
              i=i+2;
          }
          if(rep==0 && h!=myNetwAddr){
              qi[h]= 0;
              ps[h]=0;
              lamda[h]=0;
              vh[h]=0;
              r[h]=0;
              if(myNetwAddr!=9)
                uhi[myNetwAddr][h]=0;
          }
        }
      for(h=0;h<V;h++){
            for(i=0;i<N;i++){
                if(h!=myNetwAddr)
                    uhi[h][i]=0;
            }
        }

    }
}



DnLm1::~DnLm1()
{
    cancelAndDelete(0);
}

void DnLm1::handleLowerMsg(cMessage* msg)
{
    EV<<"LLLLLLLLLLOOOOOOOOOWWWWWWWWWWWWEEEEEEEEEEEERRRRRRRRRRR : "<<myNetwAddr<< endl;
    DnLm1Pkt*           netwMsg        = check_and_cast<DnLm1Pkt*>(msg);
    const LAddress::L3Type& initialSrcAddr  = netwMsg->getInitialSrcAddr();
    int repp=0;
    for(i=0;i<12;i++){
       if(myNetwAddr== netwMsg->getFinalDestAddr(i)){
           i=12;
           repp=1;
       }
    }
    if(repp==0){
        delete netwMsg;
    }else{
        const cObject* pCtrlInfo = NULL;
        DnLm1Pkt* msgCopy;
        msgCopy = netwMsg;
        if (netwMsg->getKind()==RouteEstablishment ){
            lamda[initialSrcAddr]=netwMsg->getMsgLamda();
            if(repV==0){
                qiV[initialSrcAddr]=  netwMsg->getMsgQi();
            }else{
                if(repV==1)
                    qiVV[initialSrcAddr]=netwMsg->getMsgQi();
            }

            delete msgCopy;
        }else{
            if (netwMsg->getKind()==RouteEstablishment2 ){
                i=0;
                l=0;
                while (i < 16){
                    if(tablien[i] == initialSrcAddr){
                        l=tablien[i+1];
                        i=17;
                    }else
                        i=i+2;
                }
                for (k=0; k<link; k++){
                        if(initialSrcAddr!=9)
                             xhlVoisin[initialSrcAddr][k]= netwMsg->getMsgXhl(k);
                        if (myNetwAddr==9){
                             xhl[k]=netwMsg->getMsgXhl(k);
                             xhlVoisin[initialSrcAddr][k]= netwMsg->getMsgXhl(k);
                        }
                }
                delete msgCopy;
            }else{
                if (netwMsg->getKind() == DATA && netwMsg->getFinalDestAddr(i)==9 && myNetwAddr==9) {
                    sendUp(decapsMsg(msgCopy));
                    nbDataPacketsReceived++;
                }else
                    delete msgCopy;
            }
        }

        if (pCtrlInfo != NULL)
            delete pCtrlInfo;
    }
}

void DnLm1::handleLowerControl(cMessage *msg)
{
    delete msg;
}

void DnLm1::handleUpperMsg(cMessage* msg)
{
    EV<<"UUUUUUUUUUUUUUUUUUUUPPPPPPPPPPPPPPPPEEEERRRRRRRRRRRR : "<<myNetwAddr<<endl;
    DnLm1Pkt*    routeEst   = new DnLm1Pkt(msg->getName(), RouteEstablishment);
    for(i=0;i<12;i++)
        routeEst->setFinalDestAddr(i,300);
    if (route==0){
        j=0; i=0; degre=0;
        while(tabVoisin[i]!=-1){
           k=tabVoisin[i];
           routeEst->setFinalDestAddr(j,k);
           i=i+2; j++; degre++;
        }
        if (verification==1){
            if(repV==0){
               for(i=0;i<N;i++){
                   qi[i]= qiVV[i];
               }
            }else{
                if(repV==1){
                    for(i=0;i<N;i++){
                        qi[i]= qiV[i];
                    }
                }

            }
        }else
            verification=1;


/////////////////////////////////// function return uhi (i==myNetAddr)///////////////////////////////////////////////
          for( h=0;h<V;h++){
            double sumqi = 0;
            for(l=0;l<link;l++){
                if (h!=myNetwAddr)
                    sumqi=sumqi+ail[myNetwAddr][l]*xhlVoisin[h][l];
                else
                        sumqi=sumqi+ail[myNetwAddr][l]*xhl[l];
            }
            double nhi=0;
            if ( h ==myNetwAddr)
              nhi =r[myNetwAddr];
            else{
                if( myNetwAddr == 9)
                    nhi =-r[h];
                 }
            double totalUhi= (((double)w)/(sqrt(count)))*(nhi-sumqi);
               uhi[h][myNetwAddr]=uhi[h][myNetwAddr]-totalUhi;
          }

////////////////////////////////// function return vh (h=myNetAddr)//////////////////////////////////////////
                            if (myNetwAddr!=9){
                                  double x=0;
                                  double a= ((double)sigma2)/Dh;
                                  div =  log(a) /(gam*pow(ps[myNetwAddr],0.66666667));
                                  x=(((double)w)/sqrt(count))*(r[myNetwAddr]-div);
                                  vhc=vh[myNetwAddr]-x;
                                  if(vhc>0)
                                     vh[myNetwAddr]=vhc;
                                  else
                                     vh[myNetwAddr]=0;

                            }
///////////////////////////////////////// function flambda i==mynetadress/////////////////////////////////////////////
                                   double summ=0;
                                   double sump=0;
                                   for(int l=0;l<link;l++){
                                      double  sumxhl=0;
                                      for(h=0;h<V;h++){
                                          if(h!=myNetwAddr)
                                              sumxhl=sumxhl+xhlVoisin[h][l]; //yl
                                          else
                                              sumxhl=sumxhl+xhl[l];
                                      }
                                      sump=sump+(ailp[myNetwAddr][l]*cls[l]*sumxhl); //transmission power pti
                                      summ=summ+(ailm[myNetwAddr][l]*cr*sumxhl); // reception power pri
                                    }
                                    double psi=0;
                                    if (myNetwAddr==9)
                                        psi=0;
                                    else
                                        psi=ps[myNetwAddr];
                                    double total=(qi[myNetwAddr]*bi)-sump-summ-psi;
                                    double total1=(((double)w)/sqrt(count))*total;   
                                    lamdac=lamda[myNetwAddr]-total1;
                                    if(lamdac>0)
                                      lamda[myNetwAddr]=lamdac;
                                    else
                                      lamda[myNetwAddr]=0;
                                    routeEst->setMsgLamda(lamda[myNetwAddr]);
///////////////////////////////// function wl//////////////////////////////////////
                               for( l=0;l<link;l++){
                                   double sumq=0;
                                   for ( i=0;i<N;i++){
                                      sumq=sumq+ail[i][l]*qi[i];
                                   }
                                        wl[l]=wl[l]+((((double)w)/(sqrt(count)))*sumq);

                                 }
/////////////////////////////////// function return qi/////////////////////////////////////
                           double sumailwl=0;
                           for(int l=0;l<link;l++){
                               sumailwl=sumailwl+(ail[myNetwAddr][l]*wl[l]);
                           }
                           sumqiii= sumailwl-(lamda[myNetwAddr]*bi);
                           double sumqi;
                           sumqi=-sumqiii/2;
                           if(sumqi >epsilon)
                               qi[myNetwAddr]=sumqi;
                           else
                               qi[myNetwAddr]=epsilon;
                           routeEst->setMsgQi(qi[myNetwAddr]);

                           qiV[myNetwAddr]=qi[myNetwAddr];
                           qiVV[myNetwAddr]=qi[myNetwAddr]; 

//////////////////////////////////// function return psh/////////////////////////////////////////////////
                           if(myNetwAddr!=9){
                                double a= -3*lamda[myNetwAddr];
                                double lo= ((double)sigma2)/Dh;
                                double b= sqrt((pow(3*lamda[myNetwAddr],2))+(((double)((64*deltap)))* log(lo) * vh[myNetwAddr])/ gam);
                                double c =a+b;
                                sump=((double)c)/(16*deltap);
                                double sumph=pow(sump,0.6);

                                if(sumph>epsilon)
                                  ps[myNetwAddr]=sumph;
                                else
                                  {
                                    ps[myNetwAddr]=epsilon;
                                  }

                           }
////////////////////////////////// function return rh/////////////////////////////////////////////////////
                         if (myNetwAddr!=9){
                             somDb=0;
                             for(i=0;i<N;i++)
                                 somDb=somDb+RhVoisin[i];

                                   deltarh=deltar;
                                   vhrh=vh[myNetwAddr];
                                   somun=vhrh/(2*deltar);
                                   if(somun>0)
                                     r[myNetwAddr]=somun;
                                   else
                                     r[myNetwAddr]=0;
                         }

        routeEst->setInitialSrcAddr(myNetwAddr);
        routeEst->setSrcAddr(myNetwAddr);
        routeEst->setNbHops(0);
        setDownControlInfo(routeEst, LAddress::L2BROADCAST);
        assert(static_cast<cPacket*>(msg));
        routeEst->encapsulate(static_cast<cPacket*>(msg));
        sendDown(routeEst);
        route=1;
    }else{
        if(route==1){
            DnLm1Pkt*    routeEst2   = new DnLm1Pkt(msg->getName(), RouteEstablishment2);
            for(i=0;i<12;i++)
                routeEst2->setFinalDestAddr(i,300);
            i=0; j=0; degre=0;
            while(tabVoisin[i]!=-1){
              k=tabVoisin[i];
              routeEst2->setFinalDestAddr(j,k);
              j++; degre++;
              i=i+2;
            }
////////////////////////////////////// function return xhl/////////////////////////////////////////////////////////
         if (myNetwAddr!=9){
            double xht=0;
            for( l=0;l<link;l++){
              double sum=0;
                for( i=0;i<N;i++){
                    sum=sum+lamda[i]*(cls[l]*ailp[i][l] + cr*ailm[i][l])+uhi[myNetwAddr][i]*ail[i][l];
                }
              xht=-sum/(2*deltax);
              if(xht>0){ 
                 xhl[l]=xht;
              }
              else{
                 xhl[l]=0;
              }
              routeEst2->setMsgXhl(l, xhl[l]);
            }
         }

         if(repV==0)
             repV=1;
         else
             if(repV==1)
                 repV=0;
            count++;
 
            route=0;
            routeEst2->setInitialSrcAddr(myNetwAddr);
            routeEst2->setSrcAddr(myNetwAddr);
            routeEst2->setNbHops(0);
            setDownControlInfo(routeEst2, LAddress::L2BROADCAST);
            assert(static_cast<cPacket*>(msg));
            routeEst2->encapsulate(static_cast<cPacket*>(msg));
            sendDown(routeEst2);

        }else{
            DnLm1Pkt*    pkt   = new DnLm1Pkt(msg->getName(), DATA);
            pkt->setFinalDestAddr(0,9);
            pkt->setInitialSrcAddr(myNetwAddr);
            pkt->setSrcAddr(myNetwAddr);
            pkt->setNbHops(0);
            setDownControlInfo(pkt, LAddress::L2BROADCAST);
            assert(static_cast<cPacket*>(msg));
            pkt->encapsulate(static_cast<cPacket*>(msg));
            sendDown(pkt);
            nbDataPacketsSent++;
        }
    }


}

void DnLm1::finish()
{
        recordScalar("nbDataPacketsForwarded", nbDataPacketsForwarded);
        recordScalar("nbDataPacketsReceived", nbDataPacketsReceived);
        recordScalar("nbDataPacketsSent", nbDataPacketsSent);
        recordScalar("meanNbHops", (double) nbHops / (double) nbDataPacketsReceived);

}

cMessage* DnLm1::decapsMsg(DnLm1Pkt *msg)
{
    cMessage *m = msg->decapsulate();
    setUpControlInfo(m, msg->getSrcAddr());
    nbHops = nbHops + msg->getNbHops();
    delete msg;
    return m;
}
