

cplusplus {{
#include "NetwPkt_m.h"
}}

class NetwPkt;
class LAddress::L3Type extends void;

message DnLm1Pkt extends NetwPkt
{
	int              isFlood;
	LAddress::L3Type finalDestAddr[12]; 
	LAddress::L3Type initialSrcAddr; 
 	int              nbHops;
 	double		 msgQi;
 	double		 msgLamda;
 	double		 msgXhl[22];
 	
}


