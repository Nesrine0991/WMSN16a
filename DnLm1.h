
#ifndef DnLm1_h
#define DnLm1_h

#include <map>
#include <omnetpp.h>

#include "MiXiMDefs.h"
#include "BaseNetwLayer.h"
#include "SimpleAddress.h"

class SimTracer;
class DnLm1Pkt;


double sumqiii;
double logps;
double lamdaps;
double deltarh;
double vhrh;
double somun;
double sumxhl;

class MIXIM_API DnLm1 : public BaseNetwLayer
{
private:
   DnLm1(const DnLm1&);
    DnLm1& operator=(const DnLm1&);

public:
    DnLm1()
        : BaseNetwLayer()
        , headerLength(0)
        , sinkAddress()
        , useSimTracer(false)
        , floodSeqNumber(0)
        , tracer(NULL)
        , nbDataPacketsForwarded(0)
        , nbDataPacketsReceived(0)
        , nbDataPacketsSent(0)
        , finalDestAddr(0)
        , nbHops(0)
        , msgQi(0)
        , msgLamda(0)
        , msgXhl(0)
        , trace(false), stats(false), debug(false)
    {}
    /** @brief Initialization of the module and some variables*/
    virtual void initialize(int);
    virtual void finish();
    virtual ~DnLm1();
public:

    typedef matrix<double,0,1> column_vector;
    int i;
    int l;
    int k;
    int count;
    int route;
    int link;
    int N;
    int id;
    int j;
    int degre;
    double Pi;
    double psh1;
    int** ail;
    double *qi;
    double *lamda;
    double** uhi;
    double cls[22];
    double dist[22];
    double *vh;
    double *r;
    double  *ps;
    double csll[50],alpha,beta,sum1l,sum1h,sum2l,sum2h,lamdac,sumi,div,vhc,etha,mu,sumax,cr,bi,sigma2,np,gam,deltap,deltar,deltax;
    int Dh;
    int tablien[16];
    double tabVoisin[16];
    double sumqi,lamdaa,wlll,uhii,vhh,epsilon;
    double sumf;
    double w;
   double replienEtrant;
   double nbrLienS;
   double somDb;
    double RhVoisin[10];
    double xhl[22];
    double xhlVoisin[10][22];
    int V;
    int h;
    int rep;
    int repV;
    double *wl;
    int ailp[10][22];
    int ailm[10][22];
    int verification;
    double qiV[10];
    double qiVV[10];
    double lamdaV[10];
    double xhlVoisinV[10][22];

protected:
    simtime_t firstPacketGeneration;
    simtime_t lastPacketReception;


    enum messagesTypes {
        UNKNOWN=0,
        RouteEstablishment,
        RouteEstablishment2,
        ACK,
        DATA

    };


    int headerLength;
    LAddress::L3Type sinkAddress;
    bool useSimTracer;
    unsigned long floodSeqNumber;
    SimTracer *tracer;
    long nbDataPacketsForwarded;
    long nbDataPacketsReceived;
    long nbDataPacketsSent;
    long finalDestAddr;
    long nbHops;
    long msgQi;
    long msgLamda;
    long *msgXhl;


    bool trace, stats, debug;

    /** @brief Handle messages from upper layer */
    virtual void handleUpperMsg(cMessage* msg);

    /** @brief Handle messages from lower layer */
    virtual void handleLowerMsg(cMessage* msg);

    /** @brief Handle control messages from lower layer */
    virtual void handleLowerControl(cMessage* msg);

    /** @brief Decapsulate a message */
    cMessage* decapsMsg(DnLm1Pkt *msg);

};

#endif
