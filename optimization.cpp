#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include "lp_lib.h"



#define CRNODE 36        //the number of cognitive routers
#define CHANNELSIZE  20   //the global available channels in the network
#define BASICAVAILABLECHANNEL 16  //every cognitive choose this number of available channels in global available channels
#define POWERLEVEL1 10    //the optional power levels at each cognitive router
#define POWERLEVEL2 50
#define POWERLEVEL3 100
#define RANGE1 250       //the maximal transmission range responding to different power levels
#define RANGE2 373
#define RANGE3 444
#define RADIOSIZE 3     //the number of radio interface at each cognitive router
#define REQUESTNUM 8    //the number of user requests (including converge-cast, multicast and hybrid)
#define DREQUESTNUM 24  //the number of user requests with unicast
#define SPECTRALBAND 8   //the spectrum bandwidth in Mhz
#define SINRTHE   8      //SINR
#define ALPHA 4          //path loss constant
#define BETA 0.025      //antenna parameter
#define NOISE 8e-15     //noise

#define LIST_INIT_SIZE 316
#define INCREMENT 10
#define LIST_TUPLE_SIZE 50000
#define TUPLE_INCREMENT 1000
#define ADJ_INIT_SIZE 1000
#define ADJ_INCREMENT 100
#define ARR_INIT_SIZE 1000
#define ARR_INCREMENT 100
#define SAME_CHANNEL_SIZE 5

//cognitive node
typedef struct cognitivenode
{
    double pos[2];                    //position
    char channel[CHANNELSIZE];        //可用信道数，数组中每个元素0代表该信道不可用，1代表该信道可用
    int power[RADIOSIZE];                        //可调功率等级
    char radio[RADIOSIZE];            //可用无线接口数
} CR,*PCR;
//定义最短路径上链路的结构体
typedef struct shortest_path
{
    int slink[2];
    struct shortest_path * next;
} SP,*PSP;
//定义用户请求结构体
typedef struct request
{
    int source[3];
    int dest[3];
    int	band_request;
} RQ,*PRQ;
//分解后的单源单目的
typedef struct decompose_request
{
    int source;
    int dest;
} DRQ,*PDRQ;
//该节点为最短路径经过的节点，队列的一个节点
typedef struct node
{
    int id;
    struct node * next;
} NODE,*PNODE;
//每次为源目的节点对找最短路时，都会为每个节点设置一个distance值和标记前驱节点
typedef struct linknode
{
    int distance;
    PNODE pre;
} LINKNODE,*PLINKNODE;
//队列，每个元素是一个NODE类型元素
typedef struct queue
{
    PNODE front;
    PNODE rear;
} QUEUE,*PQUEUE;
//每条链路的经过最短路径数
typedef struct passdegree
{
    int x;
    int y;
    int pdg;
} PASSDEGREE,*PPASSDEGREE;
//所有链路的经过最短路径数
typedef struct passdegreelist
{
    PPASSDEGREE pd;
    int length;
    int listsize;
} PASSDEGREELIST,*PPASSDEGREELIST;
//定义一个元组
typedef struct tuple
{
    int channel;
    int power;
    int link[2];
    int radio[2];
    double id;
} TUPLE,*PTUPLE;
//所有元组的链表
typedef struct tuplelist
{
    PTUPLE tup;
    int length;
    int listsize;
} TUPLELIST,*PTUPLELIST;
//定义4DCG的邻接表
typedef struct adjarr
{
    int length;
    int listsize;
    int * adj;
} ADJARR, *PADJARR;
//最大独立集中的元组
typedef struct idependenttuple
{
    int num;
    struct idependenttuple *next;
} ITUPLE, *PITUPLE;
//定义一个最大独立集
typedef struct independentset
{
    PITUPLE itup;
    struct independentset *inext;
} IS, *PIS;
//同一锚节点的所有独立集
typedef struct independentarr
{
    int anchor;
    PIS is;
} IARR, *PIARR;
//所有极大独立集
typedef struct isetarr
{
    int length;
    int listsize;
    PIARR s;
} ISETARR, *PISETARR;
//极大独立集实现方式二:独立集集合
typedef struct mis
{
    PITUPLE ms;
} MIS,*PMIS;
//极大独立集实现方式二:独立集指针
typedef struct istag
{
    int length;
    int listsize;
    PMIS set;
} ISTAG, *PISTAG;
//极大独立集实现方式三:独立集中元组
typedef struct ctuple
{
    int num;
    int src;
    int dest;
    int power;
} CTUPLE,*PCTUPLE;
//极大独立集实现方式三:独立集集合
typedef struct sinr
{
    int listsize;
    int length;
    PCTUPLE ct;
} SINR,*PSINR;
//极大独立集实现方式三:指向独立集集合指针
typedef struct cmis
{
    PSINR cm;
} CMIS,*PCMIS;
//极大独立集实现方式三:总索引指针
typedef struct issin
{
    int length;
    int listsize;
    PCMIS arr;
} ISSIN, *PISSIN;
//带有带宽的链路
typedef struct link
{
    int src;
    int dest;
    double bandwidth;
    struct link *next;
} LINK,*PLINK;
//数据处理后的独立集
typedef struct linkset
{
    PLINK pl;
} LINKSET,*PLINKSET;
//每条链路上的独立集
typedef struct indset
{
    int num;
    double bw;
    struct indset *next;
} INDSET,*PINDSET;
//链路索引
typedef struct linkindex
{
    PINDSET lk;
    int start;
    int end;
} LINKINDEX,*PLINKINDEX;
//传入lpsolve的请求
typedef struct source_destination
{
    int source;
    int dest;
    double bandwidth;
} SANDD,*PSANDD;
//初始化网络拓扑
void Init_topology(PCR topology)
{
    srand((unsigned)time(NULL));  //生成种子

    for(int i=0; i<CRNODE; i++)
    {
        //初始化节点位置坐标
        topology[i].pos[0] = (i/(int)sqrt((double)CRNODE))*200;
        topology[i].pos[1] = (i%(int)sqrt((double)CRNODE))*200;
        //printf("%d=(%lf,%lf)\n",i, topology[i].pos[0],topology[i].pos[1]);

        //初始化可用信道
        for(int m=0; m<CHANNELSIZE; m++)
        {
            topology[i].channel[m] = 1;
        }

        for(int j=1; j<=CHANNELSIZE-BASICAVAILABLECHANNEL; j++)
        {
            int k = rand()%CHANNELSIZE;
            topology[i].channel[k] = 0;
        }

        //初始化功率等级,一个radio一个power
        for(int n=0; n<RADIOSIZE; n++)
        {
            topology[i].power[n] = POWERLEVEL3;
        }

        //初始化可用无线接口
        for(int m=0; m<RADIOSIZE; m++)
        {
            topology[i].radio[m] = 1;
        }

    }

}
//初始化随机拓扑
void Init_random_topology(PCR topology)
{
    srand((unsigned)time(NULL));  //生成种子
    double p = 0;

    for(int i=0; i<CRNODE; i++)
    {
        //初始化节点位置坐标

        topology[i].pos[0] = rand()%1000;
        topology[i].pos[1] = rand()%1000;
        if(i>0)   //如果没有这个条件，进入while循环后就出不来了
            while(p < 50)
            {
                topology[i].pos[0] = rand()%1000;
                topology[i].pos[1] = rand()%1000;

                for(int j=0; j < i; j++)
                {
                    p = sqrt(pow((topology[i].pos[0] - topology[j].pos[0]),2) + pow((topology[i].pos[1] - topology[j].pos[1]),2));
                    if(p < 50)
                        break;
                }
            }
        // printf("%d=(%lf,%lf)\n",i, topology[i].pos[0],topology[i].pos[1]);

        //初始化可用信道
        for(int m=0; m<CHANNELSIZE; m++)
        {
            topology[i].channel[m] = 1;
        }
        for(int j=1; j<=CHANNELSIZE-BASICAVAILABLECHANNEL; j++)
        {
            int k = rand()%CHANNELSIZE;
            topology[i].channel[k] = 0;
        }
        //初始化功率等级
        for(int n=0; n<RADIOSIZE; n++)
        {
            topology[i].power[n] = POWERLEVEL3;
        }
        //初始化可用无线接口
        for(int m=0; m<RADIOSIZE; m++)
        {
            topology[i].radio[m] = 1;
        }
    }
}
//初始化邻接矩阵
void Init_arcs(PCR topology,char (&parcs)[CRNODE][CRNODE])
{
    double dis = 0;
    for(int i=0; i<CRNODE; i++)
    {
        for(int j=0; j<CRNODE; j++)
            parcs[i][j]=0;
    }
    for(int i=0; i<CRNODE; i++)
    {
        for(int j=0; j<i; j++)
        {
            dis = sqrt(pow((topology[i].pos[0]-topology[j].pos[0]),2)+pow((topology[i].pos[1]-topology[j].pos[1]),2));
            if(dis <= RANGE3)
                parcs[i][j] = 1;
            //	  printf("%lf\n",dis);
        }
        for(int j=i+1; j<CRNODE; j++)
        {
            dis = sqrt(pow((topology[i].pos[0]-topology[j].pos[0]),2)+pow((topology[i].pos[1]-topology[j].pos[1]),2));
            if(dis <= RANGE3)
                parcs[i][j] = 1;
        }
    }
    /*	for(int i=0;i<CRNODE;i++)
    	{  //输出邻接矩阵
    		printf("\n");
    		for(int j=0;j<CRNODE;j++)
    		printf("%d ",parcs[i][j]);

    	}*/
}
//产生multicast用户请求
void Gen_multi_request(PRQ req)
{
    srand((unsigned)time(NULL));  //生成种子
    for(int i=0; i<REQUESTNUM; i++)
    {
        req[i].band_request = rand()%41+60; //假设速率需求在60-100Mbps之间
        for(int j=0; j<3; j++)
        {
            req[i].source[j] = 0;
            req[i].dest[j] = 0;
        }
    }
    req[1].dest[0] = 0;
    req[1].dest[1] = 1;
    req[1].dest[2] = 2;
    req[1].source[0] = 25;
    req[2].dest[0] = 3;
    req[2].dest[1] = 4;
    req[2].dest[2] = 5;
    req[2].source[0] = 28;
    req[3].dest[0] = 0;
    req[3].dest[1] = 6;
    req[3].dest[2] = 12;
    req[3].source[0] = 9;
    req[4].dest[0] = 18;
    req[4].dest[1] = 24;
    req[4].dest[2] = 30;
    req[4].source[0] = 27;
    req[5].dest[0] = 5;
    req[5].dest[1] = 11;
    req[5].dest[2] = 17;
    req[5].source[0] = 8;
    req[6].dest[0] = 23;
    req[6].dest[1] = 29;
    req[6].dest[2] = 35;
    req[6].source[0] = 26;
    req[7].dest[0] = 30;
    req[7].dest[1] = 31;
    req[7].dest[2] = 32;
    req[7].source[0] = 7;
    req[0].dest[0] = 33;
    req[0].dest[1] = 34;
    req[0].dest[2] = 35;
    req[0].source[0] = 10;

    /*	for(int i=0;i<REQUESTNUM;i++)
    	{
    		printf("%d: ",i);
    		printf(" %d ",req[i].band_request);
    		for(int j=0;j<3;j++)
    	    printf(" (%d,%d) ",req[i].source[j],req[i].dest[j]);
    		printf("\n");
    	}*/
}
//产生converge-cast用户请求
void Gen_converge_request(PRQ req)
{
    srand((unsigned)time(NULL));  //生成种子
    for(int i=0; i<REQUESTNUM; i++)
    {
        req[i].band_request = rand()%41+60; //假设速率需求在40-100Mbps之间
        for(int j=0; j<3; j++)
        {
            req[i].source[j] = 0;
            req[i].dest[j] = 0;
        }
    }
    req[1].source[0] = 0;
    req[1].source[1] = 1;
    req[1].source[2] = 2;
    req[1].dest[0] = 25;
    req[2].source[0] = 3;
    req[2].source[1] = 4;
    req[2].source[2] = 5;
    req[2].dest[0] = 28;
    req[3].source[0] = 0;
    req[3].source[1] = 6;
    req[3].source[2] = 12;
    req[3].dest[0] = 9;
    req[4].source[0] = 18;
    req[4].source[1] = 24;
    req[4].source[2] = 30;
    req[4].dest[0] = 27;
    req[5].source[0] = 5;
    req[5].source[1] = 11;
    req[5].source[2] = 17;
    req[5].dest[0] = 8;
    req[6].source[0] = 23;
    req[6].source[1] = 29;
    req[6].source[2] = 35;
    req[6].dest[0] = 26;
    req[7].source[0] = 30;
    req[7].source[1] = 31;
    req[7].source[2] = 32;
    req[7].dest[0] = 7;
    req[0].source[0] = 33;
    req[0].source[1] = 34;
    req[0].source[2] = 35;
    req[0].dest[0] = 10;

    /*	for(int i=0;i<REQUESTNUM;i++)
    	{
    		printf("%d: ",i);
    		printf(" %d ",req[i].band_request);
    		for(int j=0;j<3;j++)
    	    printf(" (%d,%d) ",req[i].source[j],req[i].dest[j]);
    		printf("\n");
    	}*/
}
//产生hybrid用户请求
void Gen_hybrid_request(PRQ req)
{
    srand((unsigned)time(NULL));  //生成种子
    for(int i=0; i<REQUESTNUM; i++)
    {
        req[i].band_request = rand()%41+60; //假设速率需求在40-100Mbps之间
        for(int j=0; j<3; j++)
        {
            req[i].source[j] = 0;
            req[i].dest[j] = 0;
        }
    }
    req[1].source[0] = 0;
    req[1].source[1] = 1;
    req[1].source[2] = 2;
    req[1].dest[0] = 30;
    req[1].dest[1] = 31;
    req[1].dest[2] = 32;
    req[2].source[0] = 33;
    req[2].source[1] = 34;
    req[2].source[2] = 35;
    req[2].dest[0] = 3;
    req[2].dest[1] = 4;
    req[2].dest[2] = 5;
    req[3].source[0] = 0;
    req[3].source[1] = 6;
    req[3].source[2] = 12;
    req[3].dest[0] = 5;
    req[3].dest[1] = 11;
    req[3].dest[2] = 17;
    req[0].source[0] = 23;
    req[0].source[1] = 29;
    req[0].source[2] = 35;
    req[0].dest[0] = 18;
    req[0].dest[1] = 24;
    req[0].dest[2] = 30;

    /*	for(int i=0;i<REQUESTNUM;i++)
    	{
    		printf("%d: ",i);
    		printf(" %d ",req[i].band_request);
    		for(int j=0;j<3;j++)
    	    printf(" (%d,%d) ",req[i].source[j],req[i].dest[j]);
    		printf("\n");
    	}*/
}
//分解用户请求成单源单目的
void Decompose_request(PRQ req, PDRQ dreq)
{
    for(int i=0; i<REQUESTNUM; i++)
    {
        if(req[i].source[1] == 0 && req[i].source[2] == 0)
        {
            //multicast
            for(int j=0; j<3; j++)
            {
                dreq[3*i+j].source = req[i].source[0];
                dreq[3*i+j].dest = req[i].dest[j];
            }
        }
        else if(req[i].dest[1] == 0 && req[i].dest[2] == 0)
        {
            //converge-cast
            for(int j=0; j<3; j++)
            {
                dreq[3*i+j].dest = req[i].dest[0];
                dreq[3*i+j].source = req[i].source[j];
            }
        }
        else
        {
            //hybrid
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    dreq[9*i+j*3+k].source = req[i].source[j];
                    dreq[9*i+j*3+k].dest = req[i].dest[k];
                }
            }
        }

    }
    /*	for(int k=0;k<DREQUESTNUM;k++)
    	{
    		printf("(%d,%d) ",dreq[k].source,dreq[k].dest);
    	}*/
}
//初始化队列
void InitQueue(PQUEUE q)
{
    q->front = q->rear = (PNODE)malloc(sizeof(NODE));
    if(!(q->front))  exit(-1);
    q->front->next = NULL;
}
//入对列
void EnQueue(PQUEUE q,int e)
{
    PNODE p = (PNODE)malloc(sizeof(NODE));
    if(!p)  exit(-1);
    p->id = e;
    p->next = NULL;
    q->rear->next = p;
    q->rear = p;
}
//出队列
int DeQueue(PQUEUE q)
{
    PNODE p;
    if(q->front == q->rear)
        return -1;
    p = q->front->next;
    int e = p->id;
    q->front->next = p->next;
    if(q->front->next == NULL)   //保证程序健壮
        q->rear = q->front;
    free(p);
    return e;
}
//判断（i，j）链路是否已经存在最短路径链表中
char Equal(int i,int j,PSP p)
{
    while(p!=NULL)
    {
        if(i == p->slink[0] && j == p->slink[1])
            return 1;
        p = p->next;
    }
    return 0;
}
//计算最短路径，并标记经过的链路
void Caculate_shortest_path(char (&parcs)[CRNODE][CRNODE],PDRQ req,PSP * sp)
{
    LINKNODE ln[CRNODE];
    QUEUE q;
    PNODE p;
    int e;

    for(int m=0; m<DREQUESTNUM; m++)
    {
        sp[m] = NULL;   //初始化最短路径数组sp
    }

    for(int j=0; j<CRNODE; j++) //初始化节点distance和前驱
    {
        ln[j].distance = 0;
        ln[j].pre = NULL;
    }

    for(int i=0; i<DREQUESTNUM; i++)
    {
        InitQueue(&q);             //初始化q队列
        for(int j=0; j<CRNODE; j++) //初始化ln数组，每组源目的节点对会产生一个ln数组
        {
            while(ln[j].pre)       //NULL就是0
            {
                p = ln[j].pre->next;  //P指向当前想要free元素的下一个元素
                free(ln[j].pre);
                ln[j].pre = p;
            }
            ln[j].distance = 0;
            ln[j].pre = NULL;
        }

        e = req[i].source;

        while(e!=-1)    //产生ln[CRNODE]数组
        {
            for(int k=0; k<CRNODE; k++)  //检查每一个待遍历节点的邻居节点
            {

                if(parcs[e][k] == 1 && k != req[i].source)    //检查e的邻居节点，并且不是源节点，因为源节点的distance是0，跟没有遍历过的节点相同
                {
                    if(ln[k].distance == 0)         //如果该节点没有被遍历过
                    {
                        EnQueue(&q,k);                     //将e的邻居节点先入对列
                        ln[k].distance = ln[e].distance + 1;
                        PNODE n = (PNODE)malloc(sizeof(NODE));
                        if(!n) exit(-1);
                        n->next = ln[k].pre;
                        n->id = e;
                        ln[k].pre = n;
                    }
                    else
                    {
                        if(ln[k].distance > ln[e].distance)       //如果该节点被遍历过，但是也是e的邻居节点
                        {
                            PNODE n = (PNODE)malloc(sizeof(NODE));
                            if(!n) exit(-1);
                            n->next = ln[k].pre;
                            n->id = e;
                            ln[k].pre = n;
                        }
                    }
                }

            }
            e = DeQueue(&q);   //将队列中的元素检查一遍
        }

        e = req[i].dest;
        while(ln[e].distance != 0)  //生成sp[DREQUESTNUM]，如果是源和目的非连通，目的节点的distance就是0，不执行while循环
        {
            p = ln[e].pre;
            while(p != NULL)
            {
                if(ln[p->id].distance < ln[e].distance) //防止同一层上的节点有边相连
                {
                    if(!Equal(p->id,e,sp[i]))
                    {
                        EnQueue(&q,p->id);
                        PSP r = (PSP)malloc(sizeof(SP));
                        if(!r)  exit(-1);
                        r->slink[0] = p->id;
                        r->slink[1] = e;
                        r->next = sp[i];
                        sp[i] = r;
                    }
                }
                p = p->next;
            }
            e = DeQueue(&q);
        }

        /*	for(int p=0;p<CRNODE;p++)   //输出distance
        	{
        	   printf("%d ",ln[p].distance);
        	}
        	printf("\n");

        	for(int i=0;i<CRNODE;i++)   //输出ln
        	{
        		PNODE l = ln[i].pre;
        		while(l!=NULL)
        		{
        			printf("%d ",l->id);
        			l = l->next;
        		}
        		printf("\n");
        	}
        	printf("$$\n");	*/
    }
    for(int j=0; j<CRNODE; j++) //free ln[]数组中的NODE节点
    {

        while(ln[j].pre)       //NULL就是0
        {
            p = ln[j].pre->next; //P指向当前想要free元素的下一个元素
            free(ln[j].pre);
            ln[j].pre = p;
        }

    }
    /*	int number = 0;
    	for(int i=0;i<DREQUESTNUM;i++)    //输出sp
    	{
    		PSP l = sp[i];
    		while(l!=NULL)
    		{
    			printf("[%d,%d] ",l->slink[0],l->slink[1]);
    			number++;
    			l = l->next;
    		}
    		printf("##\n");
    	}
    	printf("%d",number);*/

}
//初始化链路经过最短路径数的链表
void  Init_pdlist(PPASSDEGREELIST pdlist)
{
    pdlist->pd = (PPASSDEGREE)malloc(LIST_INIT_SIZE * sizeof(PASSDEGREE));
    if(!pdlist->pd) exit(-1);
    pdlist->length = 0;
    pdlist->listsize = LIST_INIT_SIZE;
}
//将计算完的链路插入链表
void Insert_passdegreelist(PPASSDEGREELIST pdlist,int i,int j,PSP * sp)
{
    int num = 0;
    for(int k=0; k<DREQUESTNUM; k++)
    {
        PSP p =  sp[k];
        while(p!=NULL)
        {
            if(i==p->slink[0] && j==p->slink[1])
            {
                num = num + 1;
            }
            p = p->next;
        }
    }
    if((pdlist->length) >= (pdlist->listsize))
    {
        PPASSDEGREE newbase = (PPASSDEGREE) realloc(pdlist->pd,(pdlist->listsize + INCREMENT) * sizeof(PASSDEGREE));
        if(!newbase) exit(-1);
        pdlist->pd = newbase;
        pdlist->listsize += INCREMENT;
    }
    pdlist->pd[pdlist->length].x = i;
    pdlist->pd[pdlist->length].y = j;
    pdlist->pd[pdlist->length].pdg = num;
    (pdlist->length)++;


}
//计算链路经过最短路径数
void Calculate_passdegree(PPASSDEGREELIST pdlist,PSP * sp,char (&parcs)[CRNODE][CRNODE])
{
    for(int i=0; i<CRNODE; i++)
    {
        for(int j=0; j<CRNODE; j++)
        {
            if(parcs[i][j]==1)
            {
                Insert_passdegreelist(pdlist,i,j,sp);
            }
        }
    }

    /*printf("%d",pdlist->length);
    int totalnum = 0;
    for(int i=0;i<pdlist->length;i++)
    {
       printf("(%d,%d) num:%d ",pdlist->pd[i].x,pdlist->pd[i].y,pdlist->pd[i].pdg);
       totalnum = totalnum + pdlist->pd[i].pdg;
    }
    printf("%d",totalnum);*/

}
//取得经过最短路径数最大值
int Get_pdmax(PPASSDEGREELIST pdlist)
{
    int pdmax = pdlist->pd[0].pdg;
    for(int i=1; i<pdlist->length; i++)
    {
        if(pdlist->pd[i].pdg > pdmax)
        {
            pdmax = pdlist->pd[i].pdg;
        }
    }
    return pdmax;
}

void Init_tuplelist(PTUPLELIST tplst)
{
    tplst->tup = (PTUPLE) malloc(LIST_TUPLE_SIZE * sizeof(TUPLE));
    if(!tplst->tup) exit(-1);
    tplst->length = 0;
    tplst->listsize = LIST_TUPLE_SIZE;
}
//计算每个元组的权重
double Calculate_ID(int pdg,int ch,int power,int pdmax)
{
    double weight,alpha = 0.5,beta = 0.3;
    weight = alpha * pdg / pdmax + beta  + (1 - alpha - beta) * power / POWERLEVEL3;
    // printf("%4.3lf ",id);
    return weight;
}

void Insert_tuplelist(PTUPLELIST tplst,int power,int j,int m,int n,int x,int y,int pdg,int pdmax)
{
    if((tplst->length) >= (tplst->listsize))
    {
        PTUPLE newbase = (PTUPLE)realloc(tplst->tup,(tplst->listsize+TUPLE_INCREMENT)*sizeof(TUPLE));
        if(!newbase) exit(-1);
        tplst->tup = newbase;
        tplst->listsize += TUPLE_INCREMENT;
    }
    tplst->tup[tplst->length].channel = j;
    tplst->tup[tplst->length].power = power;
    tplst->tup[tplst->length].radio[0] = m;
    tplst->tup[tplst->length].radio[1] = n;
    tplst->tup[tplst->length].link[0] = x;
    tplst->tup[tplst->length].link[1] = y;
    tplst->tup[tplst->length].id = Calculate_ID(pdg,j,power,pdmax);
    (tplst->length)++;
}
//构建元组并排序
void Sort_tuple(PPASSDEGREELIST pdlist,PTUPLELIST tplst,PCR topology,int pdmax)
{
    for(int i=0; i<pdlist->length; i++) //链路
    {
        for(int j=0; j<CHANNELSIZE; j++)   //信道
        {
            if(topology[pdlist->pd[i].x].channel[j]==1 && topology[pdlist->pd[i].y].channel[j]==1)
            {
                for(int m=0; m<RADIOSIZE; m++)    //radio
                {
                    if(topology[pdlist->pd[i].x].radio[m]==1)
                    {
                        for(int n=0; n<RADIOSIZE; n++)
                        {
                            if(topology[pdlist->pd[i].y].radio[n]==1)
                            {
                                double dis = sqrt(pow((topology[pdlist->pd[i].x].pos[0]-topology[pdlist->pd[i].y].pos[0]),2) + pow((topology[pdlist->pd[i].x].pos[1]-topology[pdlist->pd[i].y].pos[1]),2));
                                for(int k=1; k<=3; k++)
                                {
                                    switch(k)    //如果距离大于power的覆盖范围，则改元组无效
                                    {
                                    case 1:
                                        if (dis < RANGE1)
                                            Insert_tuplelist(tplst,POWERLEVEL1,j,m,n,pdlist->pd[i].x,pdlist->pd[i].y,pdlist->pd[i].pdg,pdmax);
                                        break;
                                    case 2:
                                        if (dis < RANGE2)
                                            Insert_tuplelist(tplst,POWERLEVEL2,j,m,n,pdlist->pd[i].x,pdlist->pd[i].y,pdlist->pd[i].pdg,pdmax);
                                        break;
                                    case 3:
                                        if (dis < RANGE3)
                                            Insert_tuplelist(tplst,POWERLEVEL3,j,m,n,pdlist->pd[i].x,pdlist->pd[i].y,pdlist->pd[i].pdg,pdmax);
                                        break;
                                    }

                                }
                            }
                        }
                    }

                }

            }
        }

    }
    //printf("*%d*",tplst->length);
    //对元组冒泡排序
    tuple temp;
    char flag = 0;
    for(int i=0; i<tplst->length-1; i++)
    {
        for(int j=0; j<tplst->length-1-i; j++)
        {
            if(tplst->tup[j].id < tplst->tup[j+1].id)
            {
                temp = tplst->tup[j];
                tplst->tup[j] = tplst->tup[j+1];
                tplst->tup[j+1] = temp;
                flag = 1;
            }
        }
        if(flag == 0) break;
    }
}
//计算某个4DCG顶点的邻接数组
void Insert_adjarr(int node,PADJARR &arclist,int index)
{
    if((arclist[index].length) >= (arclist[index].listsize))
    {
        int * newbase = (int *)realloc(arclist[index].adj,(arclist[index].listsize+ADJ_INCREMENT)*sizeof(int));
        if(!newbase) exit(-1);
        arclist[index].adj = newbase;
        arclist[index].listsize += ADJ_INCREMENT;
    }

    int i = arclist[index].length;
    arclist[index].adj[i] = node;
    (arclist[index].length)++;

}
//初始化4DCG的邻接表数组
void Init_arclist(PADJARR &arclist,int len,TUPLELIST tplst,PCR topology)
{
    //初始化邻接表数组
    arclist = (PADJARR)malloc(len * sizeof(ADJARR));
    if(!arclist) exit(-1);
    //初始化每个顶点的邻接表
    for(int i=0; i<len; i++)
    {
        arclist[i].adj = (int *)malloc(ADJ_INIT_SIZE * sizeof(int));
        if(!arclist[i].adj) exit(-1);
        arclist[i].length = 0;
        arclist[i].listsize = ADJ_INIT_SIZE;
    }

    double dis1 = 0,dis2 = 0,dis3 = 0,dis4 = 0,s1 = 0,s2 = 0;
    for(int i=0; i<len; i++)
    {
        for(int j=i+1; j<len; j++)
        {
            //两个元组是否公用一个radio
            if(tplst.tup[i].link[0] == tplst.tup[j].link[0])
            {
                if(tplst.tup[i].radio[0] == tplst.tup[j].radio[0])
                {
                    Insert_adjarr(j,arclist,i);
                    Insert_adjarr(i,arclist,j);
                    continue;
                }
            }
            if(tplst.tup[i].link[1] == tplst.tup[j].link[1])
            {
                if(tplst.tup[i].radio[1] == tplst.tup[j].radio[1])
                {
                    Insert_adjarr(j,arclist,i);
                    Insert_adjarr(i,arclist,j);
                    continue;
                }
            }
            if(tplst.tup[i].link[1] == tplst.tup[j].link[0])
            {
                if(tplst.tup[i].radio[1] == tplst.tup[j].radio[0])
                {
                    Insert_adjarr(j,arclist,i);
                    Insert_adjarr(i,arclist,j);
                    continue;
                }
            }
            if(tplst.tup[i].link[0] == tplst.tup[j].link[1])
            {
                if(tplst.tup[i].radio[0] == tplst.tup[j].radio[1])
                {
                    Insert_adjarr(j,arclist,i);
                    Insert_adjarr(i,arclist,j);
                    continue;
                }
            }

            //两个元组是否满足SINR阈值
            if(tplst.tup[i].channel == tplst.tup[j].channel)
            {
                dis1 = sqrt(pow((topology[tplst.tup[i].link[0]].pos[0]-topology[tplst.tup[i].link[1]].pos[0]),2) + pow((topology[tplst.tup[i].link[0]].pos[1]-topology[tplst.tup[i].link[1]].pos[1]),2));
                dis2 = sqrt(pow((topology[tplst.tup[i].link[0]].pos[0]-topology[tplst.tup[j].link[1]].pos[0]),2) + pow((topology[tplst.tup[i].link[0]].pos[1]-topology[tplst.tup[j].link[1]].pos[1]),2));
                dis3 = sqrt(pow((topology[tplst.tup[j].link[0]].pos[0]-topology[tplst.tup[i].link[1]].pos[0]),2) + pow((topology[tplst.tup[j].link[0]].pos[1]-topology[tplst.tup[i].link[1]].pos[1]),2));
                dis4 = sqrt(pow((topology[tplst.tup[j].link[0]].pos[0]-topology[tplst.tup[j].link[1]].pos[0]),2) + pow((topology[tplst.tup[j].link[0]].pos[1]-topology[tplst.tup[j].link[1]].pos[1]),2));
                //printf("[(%4.0lf,%4.0lf),(%4.0lf,%4.0lf)] ",topology[tplst.tup[i].link[1]].pos[0],topology[tplst.tup[i].link[1]].pos[1],topology[tplst.tup[j].link[0]].pos[0],topology[tplst.tup[j].link[0]].pos[1]);
                s1 = BETA*pow(dis1,-4)*tplst.tup[i].power*0.001/(BETA*pow(dis2,-4)*tplst.tup[j].power*0.001 + NOISE);
                s2 = BETA*pow(dis4,-4)*tplst.tup[j].power*0.001/(BETA*pow(dis3,-4)*tplst.tup[i].power*0.001 + NOISE);

                if(s1<SINRTHE||s2<SINRTHE)
                {
                    Insert_adjarr(j,arclist,i);
                    Insert_adjarr(i,arclist,j);
                }
            }
        }
    }
    /*for(int i=0;i<len;i++)
    {
    	printf("%d ",arclist[i].length);
    }
    for(int i=0;i<100;i++)
    {
    	printf("%d ",arclist[10000].adj[i]);
    }*/
}
//初始化极大独立集数组
void Init_isarr(PISETARR isarr)
{
    isarr->s = (PIARR)malloc(ARR_INIT_SIZE * sizeof(IARR));
    if(!isarr->s) exit(-1);
    isarr->length = 0;
    isarr->listsize = ARR_INIT_SIZE;
    for(int i=0; i<ARR_INIT_SIZE; i++)
    {
        isarr->s[i].anchor = 0;
        isarr->s[i].is = NULL;
    }
}
//判断两个元组是否有边
char Is_adj(int j,int it,PADJARR arclist)
{
    //一个顶点的邻居节点一定是有序的，所以可以折半查找
    int min = 1, max = arclist[j].length;
    int mid;
    while(min <= max)
    {
        mid = (int)((min+max)/2);
        if(it == arclist[j].adj[mid])
        {
            return 1;
        }
        else if(it > arclist[j].adj[mid])
        {
            min = mid+1;
        }
        else
        {
            max = mid-1;
        }
    }
    return 0;
}

void Append_iarr(PISETARR &isarr)
{
    PIARR newbase = (PIARR)realloc(isarr->s,(isarr->listsize + ARR_INCREMENT) * sizeof(IARR));
    if(!newbase) exit(-1);
    isarr->s = newbase;
    isarr->listsize += ARR_INCREMENT;
    for(int i=isarr->length; i<isarr->listsize; i++)
    {
        isarr->s[i].anchor = 0;
        isarr->s[i].is = NULL;
    }
}
//判断是否所有顶点都被覆盖
char Is_covered(char * &cov,int len)
{
    for(int i=len-1; i>=0; i--)
    {
        if(cov[i]==0)
        {
            return 0;
            break;
        }
    }
    return 1;
}

void Calculate_MIS(TUPLELIST &tplst,PADJARR arclist,PISETARR isarr)
{
    //初始化覆盖标记数组
    char * covered = (char *)malloc(tplst.length * sizeof(char));
    if(!covered) exit(-1);
    for(int i=0; i<tplst.length; i++)
    {
        covered[i] = 0;
    }

    int p = 0,index = 0;
    char flag = 0,insert = 0;
    //定义非邻居节点的数组
    char * nonadj = (char *)malloc(tplst.length * sizeof(char));
    if(!nonadj) exit(-1);

    while(!Is_covered(covered,tplst.length))
    {
        //计算锚顶点p的非邻居节点
        for(int k=0; k<tplst.length; k++)
        {
            nonadj[k] = 0;
        }
        for(int m=0; m<arclist[p].length; m++)
        {
            nonadj[arclist[p].adj[m]] = 1;
        }
        //当动态分配的独立集数组不够时，再次分配
        if((isarr->length) >= (isarr->listsize))
        {
            Append_iarr(isarr);
        }
        isarr->s[index].anchor = p;
        (isarr->length)++;
        covered[p] = 1;
        //遍历所有未覆盖的元组
        for(int j=p+1; j<tplst.length; j++)
        {
            //该元组与p不相邻且未被覆盖
            if(nonadj[j]==0 && covered[j]==0)
            {
                PIS pis = isarr->s[index].is;
                //遍历以p为锚顶点的所有独立集
                while(pis!=NULL)
                {
                    PITUPLE pit = pis->itup;
                    //遍历包含p的某个独立集的所有元组
                    while(pit!=NULL)
                    {
                        if(Is_adj(j,pit->num,arclist))
                        {
                            flag = 1;
                            break;
                        }
                        pit = pit->next;
                    }
                    //如果j与该独立集中的所有元组都相互独立，则将j插入该独立集，标记j为覆盖顶点，并终止对j操作
                    if(flag==0)
                    {
                        insert = 1;
                        PITUPLE newtup = (PITUPLE) malloc(sizeof(ITUPLE));
                        newtup->num = j;
                        newtup->next = pis->itup;
                        pis->itup = newtup;
                        covered[j] = 1;
                        break;
                    }
                    pis = pis->inext;
                    flag = 0;
                }
                //如果j不能插入现有的独立集，则新建一个独立集，标记j为覆盖顶点
                if(insert==0)
                {
                    PIS newis = (PIS)malloc(sizeof(IS));
                    if(!newis) exit(-1);
                    PITUPLE newtuple = (PITUPLE)malloc(sizeof(ITUPLE));
                    if(!newtuple) exit(-1);
                    newis->inext = isarr->s[index].is;
                    isarr->s[index].is = newis;
                    newis->itup = newtuple;
                    newtuple->next = NULL;
                    newtuple->num = j;
                    covered[j] = 1;
                }
            }
            // printf(" %d ",p);
        }

        //从权重最大的开始，将覆盖过的元组加入独立集，以保证构成极大独立集
        flag = 0;
        insert =0;
        for(int j=0; j<tplst.length; j++)
        {
            //j是p的非邻居节点，并且已经被覆盖
            if(nonadj[j]==0 && covered[j]==1)
            {
                PIS pis = isarr->s[index].is;
                while(pis!=NULL)
                {
                    PITUPLE pit = pis->itup;
                    while(pit!=NULL)
                    {
                        if(Is_adj(j,pit->num,arclist))
                        {
                            flag = 1;
                            break;
                        }
                        pit = pit->next;
                    }
                    if(flag==0)
                    {
                        insert = 1;
                        PITUPLE newtup = (PITUPLE)malloc(sizeof(ITUPLE));
                        newtup->num = j;
                        newtup->next = pis->itup;
                        pis->itup = newtup;
                        covered[j] = 1;
                        break;
                    }
                    pis = pis->inext;
                    flag = 0;
                }
                if(insert==0 && pis==NULL)
                {
                    PIS newis = (PIS)malloc(sizeof(IS));
                    if(!newis) exit(-1);
                    PITUPLE newtuple = (PITUPLE)malloc(sizeof(ITUPLE));
                    if(!newtuple) exit(-1);
                    newis->inext = isarr->s[index].is;
                    isarr->s[index].is = newis;
                    newis->itup = newtuple;
                    newtuple->next = NULL;
                    newtuple->num = j;
                    covered[j] = 1;
                }
            }
        }
        index++;
        flag = 0;
        insert = 0;
        //free(nonadj);
        while(covered[p]==1 && p<tplst.length)
        {
            p++;
        }
        //printf("%d end ",p);
    }

    //for(int k=0;k<1003;k++)
    //	printf("%d ",p);
}
//初始化方式二极大独立集数组
void Init_istag(PISTAG istag)
{
    istag->set = (PMIS)malloc(ARR_INIT_SIZE * sizeof(MIS));
    if(!istag->set) exit(-1);
    istag->length = 0;
    istag->listsize = ARR_INIT_SIZE;
    for(int i=0; i<ARR_INIT_SIZE; i++)
    {
        istag->set[i].ms = NULL;
    }
}
//计算极大独立集方式二
void Calculate_MIS_two(TUPLELIST &tplst,PADJARR arclist,PISTAG istag)
{
    //初始化覆盖标记数组
    char * covered = (char *)malloc(tplst.length * sizeof(char));
    if(!covered) exit(-1);
    for(int i=0; i<tplst.length; i++)
    {
        covered[i] = 0;
    }

    int p = 0,index = 0;
    char flag = 0;
    //定义非邻居节点的数组
    char * nonadj = (char *)malloc(tplst.length * sizeof(char));
    if(!nonadj) exit(-1);

    while(!Is_covered(covered,tplst.length))
    {
        //计算锚顶点p的非邻居节点
        for(int k=0; k<tplst.length; k++)
        {
            nonadj[k] = 0;
        }
        for(int m=0; m<arclist[p].length; m++)
        {
            nonadj[arclist[p].adj[m]] = 1;
        }
        //当动态分配的独立集数组不够时，再次分配
        if((istag->length) >= (istag->listsize))
        {
            PMIS newbase = (PMIS)realloc(istag->set,(istag->listsize + ARR_INCREMENT) * sizeof(MIS));
            if(!newbase) exit(-1);
            istag->set = newbase;
            istag->listsize += ARR_INCREMENT;
            for(int i=istag->length; i<istag->listsize; i++)
            {
                istag->set[i].ms = NULL;
            }
        }
        PITUPLE newtup = (PITUPLE) malloc(sizeof(ITUPLE));
        newtup->num = p;
        newtup->next = istag->set[index].ms;
        istag->set[index].ms = newtup;
        (istag->length)++;
        covered[p] = 1;
        //遍历所有未覆盖的元组
        for(int j=p+1; j<tplst.length; j++)
        {
            //该元组与p不相邻且未被覆盖
            if(nonadj[j]==0 && covered[j]==0)
            {
                PITUPLE ptup = istag->set[index].ms;
                //遍历以p为锚顶点的独立集
                while(ptup!=NULL)
                {
                    if(Is_adj(j,ptup->num,arclist))
                    {
                        flag = 1;
                        break;
                    }
                    ptup = ptup->next;
                }
                //如果j与独立集中的所有元组都相互独立，则将j插入该独立集，标记j为覆盖顶点
                if(flag==0)
                {
                    PITUPLE newtup = (PITUPLE) malloc(sizeof(ITUPLE));
                    newtup->num = j;
                    newtup->next = istag->set[index].ms;
                    istag->set[index].ms = newtup;
                    covered[j] = 1;
                }
                flag = 0;
            }
        }

        //从权重最大的开始，将覆盖过的元组加入独立集，以保证构成极大独立集
        flag = 0;
        for(int j=0; j<tplst.length; j++)
        {
            //j是p的非邻居节点，并且已经被覆盖
            if(nonadj[j]==0 && covered[j]==1)
            {
                PITUPLE ptup = istag->set[index].ms;
                //遍历以p为锚顶点的独立集
                while(ptup!=NULL)
                {
                    if(Is_adj(j,ptup->num,arclist))
                    {
                        flag = 1;
                        break;
                    }
                    ptup = ptup->next;
                }
                //如果j与独立集中的所有元组都相互独立，则将j插入该独立集，标记j为覆盖顶点
                if(flag==0)
                {
                    PITUPLE newtup = (PITUPLE) malloc(sizeof(ITUPLE));
                    newtup->num = j;
                    newtup->next = istag->set[index].ms;
                    istag->set[index].ms = newtup;
                    covered[j] = 1;
                }
                flag = 0;
            }
        }
        index++;
        flag = 0;
        while(covered[p]==1 && p<tplst.length)
        {
            p++;
        }
    }
}
//输出极大独立集到EXCEL
void Output_MIS(ISETARR &isarr)
{
    PIS q;
    PITUPLE t;
    FILE *fq;
    if((fq = fopen("MIS.csv","w+"))==NULL) printf("error");

    int count_mis = 0;
    for(int i=0; i<isarr.length; i++)
    {
        fprintf(fq,"锚节点:%d,",isarr.s[i].anchor);
        q = isarr.s[i].is;
        while(q!=NULL)
        {
            fprintf(fq,"独立集,");
            count_mis ++;
            t = q->itup;
            while(t!=NULL)
            {
                fprintf(fq,"%d,",t->num);
                t = t->next;
            }
            q = q->inext;
        }
        fprintf(fq,"\n");
    }
    fclose(fq);
    printf("%d", count_mis); //独立集个数

}
//方式二输出极大独立集到EXCEL
void Output_MIS_two(ISTAG &istag)
{
    PITUPLE t;
    FILE *fq;
    if((fq = fopen("MIS2.csv","w+"))==NULL) printf("error");

    for(int i=0; i<istag.length; i++)
    {
        t = istag.set[i].ms;
        while(t != NULL)
        {
            fprintf(fq,"%d,",t->num);
            t = t->next;
        }
        fprintf(fq,"\n");
    }
    fclose(fq);
    printf("%d", istag.length); //独立集个数
}

//初始化方式三极大独立集总索引数组
void Init_issin(PISSIN issin)
{
    issin->arr = (PCMIS)malloc(ARR_INIT_SIZE * sizeof(CMIS));
    if(!issin->arr) exit(-1);
    issin->length = 0;
    issin->listsize = ARR_INIT_SIZE;
    for(int i=0; i<ARR_INIT_SIZE; i++)
    {
        issin->arr[i].cm = (PSINR)malloc(CHANNELSIZE * sizeof(SINR));
        //初始化按信道分类的独立集数组
        for(int j=0; j<CHANNELSIZE; j++)
        {
            issin->arr[i].cm[j].ct = (PCTUPLE)malloc(SAME_CHANNEL_SIZE * sizeof(CTUPLE));
            for(int k=0; k<SAME_CHANNEL_SIZE; k++)
            {
                issin->arr[i].cm[j].ct[k].num = 0;
                issin->arr[i].cm[j].ct[k].power = 0;
                issin->arr[i].cm[j].ct[k].src = 0;
                issin->arr[i].cm[j].ct[k].dest = 0;
            }
            issin->arr[i].cm[j].length = 0;
            issin->arr[i].cm[j].listsize = SAME_CHANNEL_SIZE;
        }
    }
}

char Check_SINR(SINR &s,PCR topology)
{
    double signal = 0, noise = 0, dis = 0,SINR = 0;
    for(int i= 0; i<s.length; i++)
    {
        for(int j=0; j<s.length; j++)
        {
            if(i==j)
            {
                dis = sqrt(pow((topology[s.ct[i].src].pos[0]-topology[s.ct[i].dest].pos[0]),2) + pow((topology[s.ct[i].src].pos[1]-topology[s.ct[i].dest].pos[1]),2));
                signal = BETA * pow(dis,-4) * s.ct[i].power * 0.001;
            }
            else
            {
                dis = sqrt(pow((topology[s.ct[j].src].pos[0]-topology[s.ct[i].dest].pos[0]),2) + pow((topology[s.ct[j].src].pos[1]-topology[s.ct[i].dest].pos[1]),2));
                noise = noise + BETA * pow(dis,-4) * s.ct[j].power * 0.001;
            }
        }
        SINR = signal/(noise + NOISE);
        if(SINR < 8) return 0;
    }
    return 1;
}
//追加空间
void Append_issin(PISSIN &issin)
{
    PCMIS newbase = (PCMIS)realloc(issin->arr,(issin->listsize + ARR_INCREMENT) * sizeof(CMIS));
    if(!newbase) exit(-1);
    issin->arr = newbase;
    issin->listsize += ARR_INCREMENT;
    for(int i=issin->length; i<issin->listsize; i++)
    {
        issin->arr[i].cm = (PSINR)malloc(CHANNELSIZE * sizeof(SINR));
        //初始化按信道分类的独立集数组
        for(int j=0; j<CHANNELSIZE; j++)
        {
            issin->arr[i].cm[j].ct = (PCTUPLE)malloc(SAME_CHANNEL_SIZE * sizeof(CTUPLE));
            for(int k=0; k<SAME_CHANNEL_SIZE; k++)
            {
                issin->arr[i].cm[j].ct[k].num = 0;
                issin->arr[i].cm[j].ct[k].power = 0;
                issin->arr[i].cm[j].ct[k].src = 0;
                issin->arr[i].cm[j].ct[k].dest = 0;
            }
            issin->arr[i].cm[j].length = 0;
            issin->arr[i].cm[j].listsize = SAME_CHANNEL_SIZE;
        }
    }
}
//计算极大独立集方式三，带SINR检测
void Calculate_MIS_channel(TUPLELIST tplst,PADJARR arclist,PISSIN issin,PCR topology)
{
    //初始化覆盖标记数组
    char * covered = (char *)malloc(tplst.length * sizeof(char));
    if(!covered) exit(-1);
    for(int i=0; i<tplst.length; i++)
    {
        covered[i] = 0;
    }

    int p = 0;
    char flag = 0;
    int arrindex = 0, cmindex = 0, ctindex = 0;
    //定义非邻居节点的数组
    char * nonadj = (char *)malloc(tplst.length * sizeof(char));
    if(!nonadj) exit(-1);

    while(!Is_covered(covered,tplst.length))
    {
        //计算锚顶点p的非邻居节点
        for(int k=0; k<tplst.length; k++)
        {
            nonadj[k] = 0;
        }
        for(int m=0; m<arclist[p].length; m++)
        {
            nonadj[arclist[p].adj[m]] = 1;
        }
        //当动态分配的独立集数组不够时，再次分配
        if((issin->length) >= (issin->listsize))
        {
            Append_issin(issin);
        }
        //当p所在信道分配的内存不够时，追加一个单位内存
        arrindex = issin->length;
        cmindex = tplst.tup[p].channel;
        ctindex  = issin->arr[issin->length].cm[tplst.tup[p].channel].length;

        if(issin->arr[arrindex].cm[cmindex].length >= issin->arr[arrindex].cm[cmindex].listsize)
        {
            PCTUPLE newtup = (PCTUPLE)realloc(issin->arr[arrindex].cm[cmindex].ct,(issin->arr[arrindex].cm[cmindex].listsize + 1) * sizeof(CTUPLE));
            if(!newtup) exit(-1);
            issin->arr[arrindex].cm[cmindex].ct = newtup;
            issin->arr[arrindex].cm[cmindex].listsize += 1;
        }
        //将p插入新的独立集，注意：信道取值为0-CHANNELSIZE
        issin->arr[arrindex].cm[cmindex].ct[ctindex].num = p;
        issin->arr[arrindex].cm[cmindex].ct[ctindex].power = tplst.tup[p].power;
        issin->arr[arrindex].cm[cmindex].ct[ctindex].src = tplst.tup[p].link[0];
        issin->arr[arrindex].cm[cmindex].ct[ctindex].dest = tplst.tup[p].link[1];
        (issin->arr[issin->length].cm[tplst.tup[p].channel].length)++;
        covered[p] = 1;
        //遍历所有未覆盖的元组
        for(int j=p+1; j<tplst.length; j++)
        {
            //该元组与p不相邻且未被覆盖
            if(nonadj[j]==0 && covered[j]==0)
            {

                //遍历以与j同信道的独立集
                for(int m=0; m<issin->arr[issin->length].cm[tplst.tup[j].channel].length; m++)
                {
                    if(Is_adj(j,issin->arr[issin->length].cm[tplst.tup[j].channel].ct[m].num,arclist))
                    {
                        flag = 1;
                        break;
                    }
                }

                //如果j与该独立集中的所有元组都相互独立，则将j插入该独立集，标记j为覆盖顶点
                if(flag==0)
                {

                    //当p所在信道分配的内存不够时，追加一个单位内存
                    cmindex = tplst.tup[j].channel;
                    ctindex  = issin->arr[issin->length].cm[tplst.tup[j].channel].length;

                    if(issin->arr[arrindex].cm[cmindex].length >= issin->arr[arrindex].cm[cmindex].listsize)
                    {
                        PCTUPLE newtup = (PCTUPLE)realloc(issin->arr[arrindex].cm[cmindex].ct,(issin->arr[arrindex].cm[cmindex].listsize + 1) * sizeof(CTUPLE));
                        if(!newtup) exit(-1);
                        issin->arr[arrindex].cm[cmindex].ct = newtup;
                        issin->arr[arrindex].cm[cmindex].listsize += 1;
                    }
                    //将p插入新的独立集，注意：信道取值为0-CHANNELSIZE
                    issin->arr[arrindex].cm[cmindex].ct[ctindex].num = j;
                    issin->arr[arrindex].cm[cmindex].ct[ctindex].power = tplst.tup[j].power;
                    issin->arr[arrindex].cm[cmindex].ct[ctindex].src = tplst.tup[j].link[0];
                    issin->arr[arrindex].cm[cmindex].ct[ctindex].dest = tplst.tup[j].link[1];
                    (issin->arr[issin->length].cm[tplst.tup[j].channel].length)++;
                    if(!Check_SINR(issin->arr[arrindex].cm[cmindex],topology))
                    {
                        (issin->arr[issin->length].cm[tplst.tup[j].channel].length)--;
                    }
                    else
                    {
                        covered[j] = 1;
                    }
                }

                flag = 0;
            }
        }

        //从权重最大的开始，将覆盖过的元组加入独立集，以保证构成极大独立集
        flag = 0;
        for(int j=0; j<p; j++)
        {
            //j是p的非邻居节点，并且已经被覆盖
            if(nonadj[j]==0 && covered[j]==1)
            {
                //遍历以与j同信道的独立集
                for(int m=0; m<issin->arr[issin->length].cm[tplst.tup[j].channel].length; m++)
                {
                    if(Is_adj(j,issin->arr[issin->length].cm[tplst.tup[j].channel].ct[m].num,arclist))
                    {
                        flag = 1;
                        break;
                    }
                }

                //如果j与该独立集中的所有元组都相互独立，则将j插入该独立集，标记j为覆盖顶点
                if(flag==0)
                {

                    //当p所在信道分配的内存不够时，追加一个单位内存
                    cmindex = tplst.tup[j].channel;
                    ctindex  = issin->arr[issin->length].cm[tplst.tup[j].channel].length;

                    if(issin->arr[arrindex].cm[cmindex].length >= issin->arr[arrindex].cm[cmindex].listsize)
                    {
                        PCTUPLE newtup = (PCTUPLE)realloc(issin->arr[arrindex].cm[cmindex].ct,(issin->arr[arrindex].cm[cmindex].listsize + 1) * sizeof(CTUPLE));
                        if(!newtup) exit(-1);
                        issin->arr[arrindex].cm[cmindex].ct = newtup;
                        issin->arr[arrindex].cm[cmindex].listsize += 1;
                    }
                    //将p插入新的独立集，注意：信道取值为0-CHANNELSIZE
                    issin->arr[arrindex].cm[cmindex].ct[ctindex].num = j;
                    issin->arr[arrindex].cm[cmindex].ct[ctindex].power = tplst.tup[j].power;
                    issin->arr[arrindex].cm[cmindex].ct[ctindex].src = tplst.tup[j].link[0];
                    issin->arr[arrindex].cm[cmindex].ct[ctindex].dest = tplst.tup[j].link[1];
                    (issin->arr[issin->length].cm[tplst.tup[j].channel].length)++;
                    if(!Check_SINR(issin->arr[arrindex].cm[cmindex],topology))
                    {
                        (issin->arr[issin->length].cm[tplst.tup[j].channel].length)--;
                    }

                }

                flag = 0;
            }
        }

        flag = 0;
        (issin->length)++;

        while(covered[p]==1 && p<tplst.length)
        {
            p++;
        }
    }

    /*for(int k=0;k<100;k++)
    {
    	for(int m=0;m<CHANNELSIZE;m++)
      	{
    		printf("%d  ",issin->arr[k].cm[m].length);
        }
    	printf("\n");
    }*/
}
//初始化链路带宽
void Init_linkset(PLINKSET &ls,int len)
{
    ls = (PLINKSET) malloc(len * sizeof(LINKSET));
    if(!ls) exit(-1);
    for(int i=0; i<len; i++)
    {
        ls[i].pl = NULL;
    }
}
//计算链路带宽
void Calculate_link_bandwidth(PISSIN issin,PLINKSET ls,PCR topology)
{
    //分别计算独立集中每个元组的带宽
    double signal = 0, noise = 0, dis = 0,SINR = 0,bd = 0;
    int flag = 0;
    for(int i=0; i< issin->length; i++)
    {
        for(int j=0; j< CHANNELSIZE; j++)
        {
            for(int k= 0; k<issin->arr[i].cm[j].length; k++)
            {
                //计算带宽
                noise = 0;
                for(int t= 0; t<issin->arr[i].cm[j].length; t++)
                {
                    if(k==t)
                    {
                        dis = sqrt(pow((topology[issin->arr[i].cm[j].ct[k].src].pos[0]-topology[issin->arr[i].cm[j].ct[k].dest].pos[0]),2) + pow((topology[issin->arr[i].cm[j].ct[k].src].pos[1]-topology[issin->arr[i].cm[j].ct[k].dest].pos[1]),2));
                        signal = BETA * pow(dis,-4) * issin->arr[i].cm[j].ct[k].power * 0.001;

                    }
                    else
                    {
                        dis = sqrt(pow((topology[issin->arr[i].cm[j].ct[t].src].pos[0]-topology[issin->arr[i].cm[j].ct[k].dest].pos[0]),2) + pow((topology[issin->arr[i].cm[j].ct[t].src].pos[1]-topology[issin->arr[i].cm[j].ct[k].dest].pos[1]),2));
                        noise = noise + BETA * pow(dis,-4) * issin->arr[i].cm[j].ct[t].power * 0.001;
                    }
                }

                SINR = signal/(noise + NOISE);

                bd = SPECTRALBAND * log(1+SINR)/log(2.0);
                //合并链路带宽
                flag = 0;
                PLINK temp = ls[i].pl;
                if(temp == NULL)
                {
                    PLINK t = (PLINK)malloc(sizeof(LINK));
                    if(!t)  exit(-1);
                    t->src = issin->arr[i].cm[j].ct[k].src;
                    t->dest = issin->arr[i].cm[j].ct[k].dest;
                    t->bandwidth = bd;
                    t->next = NULL;
                    ls[i].pl = t;
                    flag = 1;
                }
                else
                {
                    while(temp != NULL)
                    {
                        if((temp->src==issin->arr[i].cm[j].ct[k].src) && (temp->dest==issin->arr[i].cm[j].ct[k].dest))
                        {
                            temp->bandwidth += bd;
                            flag = 1;
                            break;
                        }
                        temp = temp->next;
                    }
                }
                if(flag == 0)
                {
                    temp = (PLINK)malloc(sizeof(LINK));
                    if(!temp)  exit(-1);
                    temp->src = issin->arr[i].cm[j].ct[k].src;
                    temp->dest = issin->arr[i].cm[j].ct[k].dest;
                    temp->bandwidth = bd;
                    temp->next = ls[i].pl;
                    ls[i].pl = temp;
                }
                flag = 0;
            }

        }

    }
    /*  for(int n=0;n<issin->length;n++)
      {
      	PLINK plnk = ls[n].pl;
      	while(plnk != NULL)
      	{
      		printf("%f ",plnk->bandwidth);
      		plnk = plnk->next;
      	}
      	printf("\n");
      }*/

}
//初始化LINKINDEX
void Init_linkindex(PLINKINDEX &links,int num)
{
    links = (PLINKINDEX) malloc(num * sizeof(LINKINDEX));
    if(!links) exit(-1);
    for(int i=0; i<num; i++)
    {
        links[i].lk = NULL;
        links[i].start = 0;
        links[i].end = 0;
    }
}
//计算链路的序号
int Calculate_link_order(char (&parcs)[CRNODE][CRNODE],PLINK lin)
{
    int count = 0;
    for(int j=0; j<CRNODE; j++)
    {
        for(int k=0; k<CRNODE; k++)
        {
            if(parcs[j][k] == 1)
                count++;
            if(j==lin->src && k==lin->dest)
                return count;
        }
    }
    return -1;
}
//转换以独立集为索引的链表为以链路为索引的链表
void Convert_to_linkindex(PLINKINDEX links,PLINKSET ls,char (&parcs)[CRNODE][CRNODE],int misnum,int linknum)
{
    int count = 0;
    char adj[CRNODE][CRNODE];
    for(int i=0; i<CRNODE; i++)
    {
        for(int j=0; j<CRNODE; j++)
        {
            adj[i][j] = parcs[i][j];
        }
    }
    for(int i=0; i<CRNODE; i++)
    {
        for(int j=0; j<CRNODE; j++)
        {
            if(parcs[i][j] == 1)
            {
                count++;
                links[count-1].start = i;
                links[count-1].end = j;
            }
        }
    }
    for(int i=0; i<misnum; i++)
    {
        PLINK lin = ls[i].pl;
        while(lin!=NULL)
        {
            count = Calculate_link_order(adj,lin);
            if(count <= linknum)
            {
                PINDSET is = (PINDSET) malloc(sizeof(INDSET));
                if(!is) exit(-1);
                is->num = i;
                is->bw = lin->bandwidth;
                is->next = links[count-1].lk;
                links[count-1].lk = is;
            }
            lin = lin->next;
        }
    }

    /* int c =0;
      for(int m=0; m<linknum; m++)
      {
          c = 0;
          printf("(%d,%d)",links[m].start,links[m].end);
          PINDSET ps = links[m].lk;
          while(ps!=NULL)
          {
              //printf("(%d,%f)",ps->num,ps->bw);
              ps = ps->next;
              c++;
          }
          printf("*%d*\n",c);
      }*/
}
//产生converge源/目的请求
void Init_converge_sandd(PRQ req, PSANDD sd, char (&parcs)[CRNODE][CRNODE], char (&newarcs)[CRNODE+REQUESTNUM][CRNODE+REQUESTNUM] )
{
    for(int i=0; i<CRNODE+REQUESTNUM; i++)
    {
        for(int j=0; j<CRNODE+REQUESTNUM; j++)
        {
            newarcs[i][j] = 0;
        }
    }
    for(int i=0; i<CRNODE; i++)
    {
        for(int j=0; j<CRNODE; j++)
        {
            newarcs[i][j] = parcs[i][j];
        }
    }

    for(int i=CRNODE; i<CRNODE+REQUESTNUM; i++)
    {
        newarcs[i][req[i-CRNODE].source[0]] = 1;
        newarcs[i][req[i-CRNODE].source[1]] = 1;
        newarcs[i][req[i-CRNODE].source[2]] = 1;
    }
    for(int i=0; i<REQUESTNUM; i++)
    {
        sd[i].source = i+CRNODE;
        sd[i].dest = req[i].dest[0];
        sd[i].bandwidth = req[i].band_request;
    }

    /*for(int i=0;i<CRNODE+REQUESTNUM;i++)
    {
        for(int j=0;j<CRNODE+REQUESTNUM;j++)
        {
            printf("%d",newarcs[i][j]);
        }
        printf("\n");
    }
    for(int i=0;i<REQUESTNUM;i++)
    {
        printf("(%d,%d,%f)",sd[i].source,sd[i].dest,sd[i].bandwidth);
    }*/

}
//产生multicast源/目的请求
void Init_multicast_sandd(PRQ req,PSANDD sd)
{
    for(int i=0; i<REQUESTNUM; i++)
    {
        for(int j=0; j<3; j++)
        {
            sd[3*i+j].dest = req[i].dest[j];
            sd[3*i+j].source = req[i].source[0];
            sd[3*i+j].bandwidth = req[i].band_request;
        }
    }
    /*for(int i=0;i<DREQUESTNUM;i++)
    {
        printf("source:%d,dest:%d,bandwidth:%f\n",sd[i].source,sd[i].dest,sd[i].bandwidth);
    }*/

}
//lpsolve求解converge
int Solve_converge(char (&parcs)[CRNODE+REQUESTNUM][CRNODE+REQUESTNUM],PLINKINDEX links,SANDD sd[], RQ rq[], int linknum,int misnum)
{

    lprec *lp;
    int Ncol, *colno = NULL, j, ret = 0;
    int alphabase, lambdabase,totalbase;
    REAL *row = NULL;
    /*Ncol表示总列数，即变量数，colno指示每一行的列数，j表示行数，row存储每列的系数*/

    Ncol = REQUESTNUM * linknum + misnum + 1 + REQUESTNUM*4; //列数
    alphabase = REQUESTNUM*linknum;
    lambdabase = alphabase+misnum+1;
    totalbase = lambdabase+REQUESTNUM;

    lp = make_lp(0, Ncol);
    if(lp == NULL)
        ret = 1; /* couldn't construct a new model... */

    if(ret == 0)
    {
        /* create space large enough for one row */
        colno = (int *) malloc(Ncol * sizeof(*colno));
        row = (REAL *) malloc(Ncol * sizeof(*row));
        if((colno == NULL) || (row == NULL))
            ret = 2;
    }
    if(ret == 0)
    {
        set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

        j = 0;
        colno[j] = lambdabase; /* first column */
        row[j++] = 1;

        /* set the objective in lpsolve */
        if(!set_obj_fnex(lp, j, row, colno))
            ret = 3;
    }

    if(ret == 0)
    {
        //每个请求的比例系数lambda满足的约束
        for(int u=0; u<REQUESTNUM; u++)
        {
            j = 0;
            colno[j] = lambdabase; /* first column */
            row[j++] = 1;

            colno[j] = lambdabase+u+1; /* second column */
            row[j++] = -1;

            /* add the row to lpsolve */
            if(!add_constraintex(lp, 2, row, colno, LE, 0))
                ret = 4;
        }
    }

    if(ret == 0)
    {
        //每个独立集系数满足的约束
        for(int p=0; p<misnum; p++)
        {
            colno[p] = alphabase+p+1; /* first column */
            row[p] = 1;
        }
        /* add the row to lpsolve */
        if(!add_constraintex(lp, misnum, row, colno, LE, 1))
            ret = 4;
    }

    if(ret == 0)
    {
        //每条链路满足的约束
        for(int i=0; i<linknum; i++)
        {
            for(j=0; j<REQUESTNUM; j++)
            {
                colno[j] = i*REQUESTNUM+j+1;
                row[j] = 1;
            }
            PINDSET li = links[i].lk;
            while(li!=NULL)
            {
                colno[j] = alphabase+(li->num)+1;
                row[j] = -li->bw;
                li = li->next;
                j++;
            }
            if(!add_constraintex(lp, j, row, colno, LE, 0))
                ret = 4;
        }
    }

    if(ret == 0)
    {

        for(int i=0; i<REQUESTNUM; i++)
        {
            //源节点发出的流等于用户请求
            j = 0;

            colno[j] = totalbase+i*3+1;
            row[j++] = 1;

            colno[j] = totalbase+i*3+2;
            row[j++] = 1;

            colno[j] = totalbase+i*3+3;
            row[j++] = 1;

            colno[j] = lambdabase+i+1;
            row[j++] = -sd[i].bandwidth;
            if(!add_constraintex(lp, j, row, colno, EQ, 0))
                ret = 4;
            //进入源节点的流为0
            j = 0;
            int t = 0;

            for(int s=0; s<3; s++)
            {
                t = rq[i].source[s];
                for(int k=0; k<CRNODE; k++)
                {
                    if(parcs[k][t]==1)
                    {
                        int s = 0;
                        for(int m=0; m<linknum; m++)
                        {
                            if(k==links[m].start && t==links[m].end)
                            {
                                s = m;
                                break;
                            }
                        }
                        colno[j] = s*REQUESTNUM+i+1;
                        row[j++] = 1;
                    }
                }

            }
            if(!add_constraintex(lp, j, row, colno, EQ, 0))
                ret = 4;
            //目的节点发出的流为0
            j = 0;
            t = sd[i].dest;
            for(int k=0; k<CRNODE; k++)
            {
                if(parcs[t][k]==1)
                {
                    int s = 0;
                    for(int m=0; m<linknum; m++)
                    {
                        if(t==links[m].start && k==links[m].end)
                        {
                            s = m;
                            break;
                        }
                    }
                    colno[j] = s*REQUESTNUM+i+1;
                    row[j++] = 1;
                }
            }
            if(!add_constraintex(lp, j, row, colno, EQ, 0))
                ret = 4;
            //中间节点流入等于流出
            j = 0;
            for(int k=0; k<CRNODE; k++)
            {
                if(k!=rq[i].source[0]&&k!=rq[i].source[1]&&k!=rq[i].source[2] && k!=sd[i].dest)
                {
                    for(int p=0; p<CRNODE; p++)
                    {
                        if(parcs[p][k]==1)
                        {
                            int s = 0;
                            for(int m=0; m<linknum; m++)
                            {
                                if(p==links[m].start && k==links[m].end)
                                {
                                    s = m;
                                    break;
                                }
                            }
                            colno[j] = s*REQUESTNUM+i+1;
                            row[j++] = 1;
                        }
                    }

                    for(int q=0; q<CRNODE; q++)
                    {
                        if(parcs[k][q]==1)
                        {
                            int s = 0;
                            for(int m=0; m<linknum; m++)
                            {
                                if(k==links[m].start && q==links[m].end)
                                {
                                    s = m;
                                    break;
                                }
                            }
                            colno[j] = s*REQUESTNUM+i+1;
                            row[j++] = -1;
                        }
                    }
                    if(!add_constraintex(lp, j, row, colno, EQ, 0))
                        ret = 4;
                    j = 0;
                }
            }
            //真实源节点的处理
            j = 0;
            for(int r=0; r<3; r++)
            {
                t = rq[i].source[r];
                for(int q=0; q<CRNODE; q++)
                {
                    if(parcs[t][q]==1)
                    {
                        int s = 0;
                        for(int m=0; m<linknum; m++)
                        {
                            if(t==links[m].start && q==links[m].end)
                            {
                                s = m;
                                break;
                            }
                        }
                        colno[j] = s*REQUESTNUM+i+1;
                        row[j++] = 1;
                    }
                }
                colno[j] = totalbase+i*3+r;
                row[j++] = -1;
                if(!add_constraintex(lp, j, row, colno, EQ, 0))
                    ret = 4;
            }
        }

    }
    if(ret == 0)
    {
        //每个lambda<=1
        for(int u=0; u<REQUESTNUM; u++)
        {
            j = 0;
            colno[j] = lambdabase+u+1; /* first column */
            row[j++] = 1;
            /* add the row to lpsolve */
            if(!add_constraintex(lp, j, row, colno, LE, 1))
                ret = 4;
        }
        set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
    }

    if(ret == 0)
    {
        /* set the object direction to maximize */
        set_maxim(lp);

        /* just out of curioucity, now show the model in lp format on screen */
        /* this only works if this is a console application. If not, use write_lp and a filename */
        // write_LP(lp, stdout);
        write_lp(lp, "model.lp");

        /* I only want to see important messages on screen while solving */
        set_verbose(lp, IMPORTANT);

        /* Now let lpsolve calculate a solution */
        ret = solve(lp);
        if(ret == OPTIMAL)
            ret = 0;
        else
            ret = 5;
    }

    if(ret == 0)
    {
        /* a solution is calculated, now lets get some results */

        /* objective value */
        printf("Objective value: %f\n", get_objective(lp));

        /* variable values */
        /* get_variables(lp, row);
         for(j = 0; j < Ncol; j++)
           printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);*/

        /* we are done now */
    }

    /* free allocated memory */
    if(row != NULL)
        free(row);
    if(colno != NULL)
        free(colno);

    if(lp != NULL)
    {
        /* clean up such that all used memory by lpsolve is freed */
        delete_lp(lp);
    }

    return(ret);
}
//lpsolve求解multicast
int Solve_multicast(char (&parcs)[CRNODE][CRNODE],PLINKINDEX links,SANDD rq[],int linknum,int misnum)
{

    lprec *lp;
    int Ncol, *colno = NULL, j, ret = 0;
    int alphabase, lambdabase;
    REAL *row = NULL;
    /*Ncol表示总列数，即变量数，colno指示每一行的列数，j表示行数，row存储每列的系数*/

    Ncol = DREQUESTNUM * linknum + misnum + 1 + DREQUESTNUM; //列数
    alphabase = DREQUESTNUM*linknum;
    lambdabase = alphabase+misnum+1;

    lp = make_lp(0, Ncol);
    if(lp == NULL)
        ret = 1; /* couldn't construct a new model... */

    if(ret == 0)
    {
        /* create space large enough for one row */
        colno = (int *) malloc(Ncol * sizeof(*colno));
        row = (REAL *) malloc(Ncol * sizeof(*row));
        if((colno == NULL) || (row == NULL))
            ret = 2;
    }
    if(ret == 0)
    {
        set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

        j = 0;
        colno[j] = lambdabase; /* first column */
        row[j++] = 1;

        /* set the objective in lpsolve */
        if(!set_obj_fnex(lp, j, row, colno))
            ret = 3;
    }

    if(ret == 0)
    {

        for(int u=0; u<DREQUESTNUM; u++)
        {
            j = 0;
            colno[j] = lambdabase; /* first column */
            row[j++] = 1;

            colno[j] = lambdabase+u+1; /* second column */
            row[j++] = -1;

            /* add the row to lpsolve */
            if(!add_constraintex(lp, 2, row, colno, LE, 0))
                ret = 4;
        }

    }

    if(ret == 0)
    {

        for(int p=0; p<misnum; p++)
        {
            colno[p] = alphabase+p+1; /* first column */
            row[p] = 1;
        }
        /* add the row to lpsolve */
        if(!add_constraintex(lp, misnum, row, colno, LE, 1))
            ret = 4;
    }

    if(ret == 0)
    {

        for(int i=0; i<linknum; i++)
        {
            for(j=0; j<DREQUESTNUM; j++)
            {
                colno[j] = i*DREQUESTNUM+j+1;
                row[j] = 1;
            }
            PINDSET li = links[i].lk;
            while(li!=NULL)
            {
                colno[j] = alphabase+(li->num)+1;
                row[j] = -li->bw;
                li = li->next;
                j++;
            }
            if(!add_constraintex(lp, j, row, colno, LE, 0))
                ret = 4;
        }
    }

    if(ret == 0)
    {

        for(int i=0; i<DREQUESTNUM; i++)
        {
            //源节点发出的流等于用户请求
            j = 0;
            int t = rq[i].source;
            for(int k=0; k<CRNODE; k++)
            {
                if(parcs[t][k]==1)
                {
                    int s = 0;
                    for(int m=0; m<linknum; m++)
                    {
                        if(t==links[m].start && k==links[m].end)
                        {
                            s = m;
                            break;
                        }
                    }
                    colno[j] = s*DREQUESTNUM+i+1;
                    row[j++] = 1;
                }
            }
            colno[j] = lambdabase+i+1;
            row[j++] = -rq[i].bandwidth;
            if(!add_constraintex(lp, j, row, colno, EQ, 0))
                ret = 4;
            //进入源节点的流为0
            j = 0;
            for(int k=0; k<CRNODE; k++)
            {
                if(parcs[k][t]==1)
                {
                    int s = 0;
                    for(int m=0; m<linknum; m++)
                    {
                        if(k==links[m].start && t==links[m].end)
                        {
                            s = m;
                            break;
                        }
                    }
                    colno[j] = s*DREQUESTNUM+i+1;
                    row[j++] = 1;
                }
            }
            if(!add_constraintex(lp, j, row, colno, EQ, 0))
                ret = 4;
            //目的节点发出的流为0
            j = 0;
            t = rq[i].dest;
            for(int k=0; k<CRNODE; k++)
            {
                if(parcs[t][k]==1)
                {
                    int s = 0;
                    for(int m=0; m<linknum; m++)
                    {
                        if(t==links[m].start && k==links[m].end)
                        {
                            s = m;
                            break;
                        }
                    }
                    colno[j] = s*DREQUESTNUM+i+1;
                    row[j++] = 1;
                }
            }
            if(!add_constraintex(lp, j, row, colno, EQ, 0))
                ret = 4;
            //中间节点流入等于流出
            j = 0;
            for(int k=0; k<CRNODE; k++)
            {
                if(k!=rq[i].source && k!=rq[i].dest)
                {
                    for(int p=0; p<CRNODE; p++)
                    {
                        if(parcs[p][k]==1)
                        {
                            int s = 0;
                            for(int m=0; m<linknum; m++)
                            {
                                if(p==links[m].start && k==links[m].end)
                                {
                                    s = m;
                                    break;
                                }
                            }
                            colno[j] = s*DREQUESTNUM+i+1;
                            row[j++] = 1;
                        }
                    }

                    for(int q=0; q<CRNODE; q++)
                    {
                        if(parcs[k][q]==1)
                        {
                            int s = 0;
                            for(int m=0; m<linknum; m++)
                            {
                                if(k==links[m].start && q==links[m].end)
                                {
                                    s = m;
                                    break;
                                }
                            }
                            colno[j] = s*DREQUESTNUM+i+1;
                            row[j++] = -1;
                        }
                    }
                    if(!add_constraintex(lp, j, row, colno, EQ, 0))
                        ret = 4;
                        j = 0;
                }
            }

        }



    }
        if(ret == 0)
    {

        for(int u=0; u<DREQUESTNUM; u++)
        {
            j = 0;
            colno[j] = lambdabase+u+1; /* first column */
            row[j++] = 1;


            /* add the row to lpsolve */
            if(!add_constraintex(lp, j, row, colno, LE, 1))
                ret = 4;
        }
          set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
    }

    if(ret == 0)
    {
        /* set the object direction to maximize */
        set_maxim(lp);

        /* just out of curioucity, now show the model in lp format on screen */
        /* this only works if this is a console application. If not, use write_lp and a filename */
        // write_LP(lp, stdout);
        write_lp(lp, "model.lp");

        /* I only want to see important messages on screen while solving */
        set_verbose(lp, IMPORTANT);

        /* Now let lpsolve calculate a solution */
        ret = solve(lp);
        if(ret == OPTIMAL)
            ret = 0;
        else
            ret = 5;
    }

    if(ret == 0)
    {
        /* a solution is calculated, now lets get some results */

        /* objective value */
        printf("Objective value: %f\n", get_objective(lp));

        /* variable values */
        /* get_variables(lp, row);
         for(j = 0; j < Ncol; j++)
           printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);*/

        /* we are done now */
    }

    /* free allocated memory */
    if(row != NULL)
        free(row);
    if(colno != NULL)
        free(colno);

    if(lp != NULL)
    {
        /* clean up such that all used memory by lpsolve is freed */
        delete_lp(lp);
    }

    return(ret);
}

int main(void)
{

    //初始化网络拓扑
    CR topology[CRNODE];
//	Init_random_topology(topology);    //随机拓扑
    Init_topology(topology);           //网格拓扑
    //初始化邻接矩阵
    char arcs[CRNODE][CRNODE];
    Init_arcs(topology,arcs);
    //count为网络中的链路数，包含双向，所以一定是偶数
    int count = 0;
    for(int i=0; i<CRNODE; i++)
    {
        for(int j=0; j<CRNODE; j++)
            if(arcs[i][j]==1) count++;
    }

    //初始化用户需求
    RQ req[REQUESTNUM];
   // Gen_multi_request(req);          //multicast
    Gen_converge_request(req);         //converge-cast
     //Gen_hybrid_request(req);        //hybrid
    //分解用户请求为单播
    DRQ dreq[DREQUESTNUM];
    Decompose_request(req,dreq);
    //计算源和目的结点间最短路径
    PSP sp[DREQUESTNUM];
    Caculate_shortest_path(arcs,dreq,sp);
    //计算每条链路的经过最短路径值
    PASSDEGREELIST pdlist;
    Init_pdlist(&pdlist);
    Calculate_passdegree(&pdlist,sp,arcs);
    int pdmax = Get_pdmax(&pdlist);
    //printf("%d",pdmax);

    //对元组排序
    TUPLELIST tplst;
    Init_tuplelist(&tplst);
    Sort_tuple(&pdlist,&tplst,topology,pdmax);
    free(pdlist.pd);

    //输出排好序的元组到文件
    /*FILE *fp;
    if((fp = fopen("sorted_tuple.csv","w+")) == NULL) printf("error");
    for(int i=0;i<tplst.length;i++)
    {
    	fprintf(fp,"%d,%d,%d,%d,%d,%d,%lf\n",tplst.tup[i].channel,tplst.tup[i].power,tplst.tup[i].link[0],tplst.tup[i].link[1],tplst.tup[i].radio[0],tplst.tup[i].radio[1],tplst.tup[i].id);
    }
    fclose(fp);*/
    //建立4DCG邻接表
    PADJARR arclist = NULL;
    Init_arclist(arclist,tplst.length,tplst,topology);

    printf("*%d*",count);

    //极大独立集实现方式三
    ISSIN issin;
    Init_issin(&issin);
    Calculate_MIS_channel(tplst,arclist,&issin,topology);
    // printf("#%d#\n ",issin.length);


    //计算每个独立集中各条链路的带宽
    PLINKSET ls = NULL;
    Init_linkset(ls,issin.length);
    Calculate_link_bandwidth(&issin,ls,topology);
    //释放之前主动分配的内存
    free(tplst.tup);
    for(int i=0; i<tplst.length; i++)
    {
        free(arclist[i].adj);
    }
    free(arclist);
    int MISnum = issin.length;
    for(int j=0; j<MISnum; j++)
    {
        for(int k=0; k<CHANNELSIZE; k++)
        {
            free(issin.arr[j].cm[k].ct);
        }
        free(issin.arr[j].cm);
    }
    free(issin.arr);
    //将LINKSET转化为以链路为索引的链表
    PLINKINDEX links = NULL;
    Init_linkindex(links,count);
    Convert_to_linkindex(links,ls,arcs,MISnum,count);

    //converge lpsolve
    SANDD sd[REQUESTNUM];
    char ArcsConverge[CRNODE+REQUESTNUM][CRNODE+REQUESTNUM];
    Init_converge_sandd(req,sd,arcs,ArcsConverge);
    Solve_converge(ArcsConverge,links,sd,req,count,MISnum);
    //multicast lpsolve
   /* SANDD sd[DREQUESTNUM];
    Init_multicast_sandd(req,sd);
    Solve_multicast(arcs,links,sd,count,MISnum);*/
    //hybrid lpsolve


    system("pause");
    return 0;

}
