typedef struct ContourNodeStruct
{
	int number;
	ContourNodeStruct* next;
	
	int x;
	int y;
	double curvature;
	double rotation;
	int distance;
	int peak=0;
}ContourNode;

typedef struct MatchingStruct
{
	int Sample_Number=0;
	CvPoint* Sample_Point=NULL;
	CvPoint Shift_Point;		//=center
/**Below are data recording angle and trigonometry **/
	double angle=0;
	double cs=1;
	double sn=0;
	double times=1;
	int min_dist;		//Initialized as infinity large
	
	MatchingStruct* Template=NULL;	//If Template=NULL, that means this seq is a template. Otherwise, this point can refer to the template.
	float AxisParam[4]={0,0,0,0};	//two vectors which indicate the center and the direction of the contour axis.
	CvPoint Big[2];					//two points which indicate the widest points of the contour
	CvPoint Fin_Left;
	CvPoint Fin_Right;
	
	
	
}Data_Match;

class SuperContour
{
	public:
/*******************************
数据存储区
*******************************/		
		CvMemStorage* Seq_Storage=NULL;
		CvSeq* OriginContour=NULL;			//处理的轮廓序列
		ContourNode* Contourlist=NULL;		//生成的轮廓序列表
		int length=0;						//轮廓长度
		CvPoint center;
		Data_Match* Match_Data=NULL;
/*******************************
数据列表操作基本函数
*******************************/
	
		void Super_Init();
		
		
		void List_Create(CvSeq* OriginContour);
		void List_Insert(int place, ContourNode* NewNode);	//If NewNode== NULL, then create an empty node with data zeroed.
		void List_Del(int place);
		void List_Replace(int place, ContourNode* NewNode);
		void List_Visit(int place, ContourNode* output);
		ContourNode* List_Use(int place);
		void List_Clear();
/*******************************
数据计算和分析函数
*******************************/
		void Contour_curvature(int N);
		void Contour_rotation(int N);
		CvPoint Contour_Center();
		void Contour_Big();
		void Contour_Distance();
/*******************************
Below are functions for matching two seqs.
*******************************/
		void Match_Sampling(int SpNumber);
		void Match_Angle(SuperContour Contour_Target,int SpNumber, double Rotate_angle);
		void Match_Best(SuperContour Contour_Target,int SpNumber);
		void Match_FitAxis();
		void Fish_NormParam(SuperContour *Template);
		
		void Transform(float cs, float sn, CvPoint shift, double times);
		void Normalize();
		void Match_Fin();
		//void Match_Rotate(double);
};

/***************************************************************/

void SuperContour::Super_Init()
{
	Seq_Storage=cvCreateMemStorage(0);
	Match_Data=new Data_Match;
	Match_Data->min_dist=1000000;
	Match_Data->Fin_Left=cvPoint(0,0);
	Match_Data->Fin_Right=cvPoint(0,0);
	
	OriginContour=NULL;
	Contourlist=NULL;
	length=0;
	center=cvPoint(0,0);
}

void SuperContour::List_Create(CvSeq* Origin)
{
	length=0;
	ContourNode* newnode=NULL;
	CvPoint* p=NULL;
	
	if(Contourlist!=NULL)
	{
		printf("\n The LIST is not empty.");
		return;
	}
	if(Origin==NULL)
	{
		printf("The ORIGIN is empty, can not copy data in.");
		return;
	}
	OriginContour=Origin;
	for(int i=0;i<Origin->total;i++)
	{
		///printf("YES\n");
		p=(CvPoint*)cvGetSeqElem(Origin,i);
		if(i==0)
		{
			newnode=new ContourNode;
			Contourlist=newnode;
		}		
		else
		{
			newnode->next=new ContourNode;
			newnode=newnode->next;
		}
		
		newnode->number=i;
		length++;
		
		newnode->x=p->x;
		newnode->y=p->y;
		///printf("\nReading point no.(%d), (%d,%d)\n",i,p->x,p->y);
		newnode->curvature=0.0;
		newnode->rotation=0.0;
		
	}
	newnode->next=NULL;
	center=Contour_Center();
	
}

void SuperContour::List_Insert(int place, ContourNode* NewNode)
{
	if(NewNode==NULL)
	{
		printf("Empty NEWNODE, return.");
		return;
	}
	ContourNode *newnode=NULL,*tempnode=NULL;
	if(Contourlist==NULL)
	{
		printf("\nThe LIST is empty, create 1 node now.");
		newnode=new ContourNode;
	
		newnode->number=0;
		newnode->x=NewNode->x;
		newnode->y=NewNode->y;
		
		newnode->curvature=0;
		
		newnode->next=NULL;
		Contourlist=newnode;
	}
	else
	{
		ContourNode* nextnode;
		for(nextnode=Contourlist->next;nextnode->next!=NULL;nextnode=nextnode->next);
		{
		}
		newnode=new ContourNode;
	
		newnode->number=0;
		newnode->x=NewNode->x;
		newnode->y=NewNode->y;
		
		newnode->curvature=0;
		
		newnode->next=NULL;
		nextnode->next=newnode;
	}
}

void SuperContour::List_Del(int place)
{
	ContourNode *nextnode=NULL,*tempnode=NULL;
	if(Contourlist==NULL)
	{
		printf("\nEmpty LIST.");
		return;
	}
	
	if(place<0)
	{
		printf("\nInvalid input INDEX.");
		return;
	}
	
	length--;
	
	if(place==0)
		{			
			tempnode=Contourlist;
			Contourlist=Contourlist->next;
			delete tempnode;
			tempnode=NULL;
			return;
		}
	for(nextnode=Contourlist;nextnode->next!=NULL;)
	{		
		if(place==nextnode->number)
		{
			tempnode->next=nextnode->next;
			delete nextnode;
			nextnode=NULL;
			break;
		}
		tempnode=nextnode;
		nextnode=nextnode->next;
	}
	if(tempnode->next==NULL)
	{
		return;
	}
	for(;tempnode!=NULL;tempnode=tempnode->next)
	{
		tempnode->number--;
	}	
	
}

void SuperContour::List_Replace(int place, ContourNode* NewNode)
{
	ContourNode *tempnode=NULL;
	if(Contourlist==NULL||NewNode==NULL)
	{
		printf("Empty INPUT (LIST/Node).");
		return;
	}
	
	for(tempnode=Contourlist;tempnode!=NULL;tempnode=tempnode->next)
	{
		if(tempnode->number==place)
		{
			tempnode->x=NewNode->x;
			tempnode->y=NewNode->y;
			
			tempnode->curvature;
			return;
		}
	}
	
}

void SuperContour::List_Visit(int place, ContourNode* output)
{
	ContourNode* temp;
	if(place>=length||place<0)
	{
		printf("Index exceeds.\n");
		return;
	}
	if(Contourlist==NULL)
	{
		printf("Vanished Contour.\n");
		return;
	}
	if(output==NULL)
	{
		printf("Empty Pointer for data-storage.\n");
		return;
	}
	
	for(temp=Contourlist;temp!=NULL;temp=temp->next)
	{
		if(temp->number==place)
		{
			///printf("\nThe index is %d, the temp->number is %d",place,temp->number);
			
			output->number=temp->number;
			output->next=NULL;
			
			output->x=temp->x;
			output->y=temp->y;
			
			
			output->curvature=temp->curvature;
			output->rotation=temp->rotation;
			
			///printf("   output number:%d -curve: %f\n",temp->number,output->curvature);
			return;
		}
	}
	
	printf("Sth's wrong");
}

ContourNode* SuperContour::List_Use(int place)
{
	ContourNode* temp;
	ContourNode* output;
	if(place>=length||place<0)
	{
		printf("Index exceeds.\n");
		return NULL;
	}
	if(Contourlist==NULL)
	{
		printf("Vanished Contour.\n");
		return NULL;
	}
	
	
	for(temp=Contourlist;temp!=NULL;temp=temp->next)
	{
		if(temp->number==place)
		{
			///printf("\nThe index is %d, the temp->number is %d",place,temp->number);
			
			
			output=temp;
			///printf("   output number:%d -curve: %f\n",temp->number,output->curvature);
			return output;
		}
	}
	
	printf("Sth's wrong");
}

void SuperContour::List_Clear()
{
	ContourNode* temp1;
	ContourNode* temp2;
	
	if(Contourlist==NULL)
	{
		length=0;
		delete OriginContour;
		OriginContour=NULL;
		printf("All clear.\n");
		return;
	}
	if(Contourlist->next==NULL)
	{
		length=0;
		delete Contourlist;
		Contourlist=NULL;
		delete OriginContour;
		OriginContour=NULL;
		printf("All clear.\n");
		return;
	}
	
	for(temp1=Contourlist;temp1->next!=NULL;)
	{
		temp2=temp1;
		temp1=temp1->next;
		delete temp2;
		temp2=NULL;
	}
	length=0;
	delete OriginContour;
	OriginContour=NULL;
	Contourlist=NULL;
}
/***************************************************************/

void SuperContour::Contour_curvature(int N)
{
	if(N>=length)
	{
		printf("\nThe LENGTH of the contour is too small, won't calculate the curvature.\n");
		return;
	}
	if(Contourlist==NULL)
	{
		printf("Empty LIST.\n");
		return;
	}
	double sum_x,sum_y,sum_xx,sum_yy,sum_xy,sum_xxx,sum_yyy,sum_xxy,sum_xyy;
	double C,D,E,G,H,a,b,c,R;
	
	int *x=(int*)malloc(length*sizeof(int));
	int *y=(int*)malloc(length*sizeof(int));
	
	ContourNode* temp=new ContourNode;
	
	for(int i=0;i<length;i++)
	{
		List_Visit(i,temp);
		x[i]=temp->x;
		y[i]=temp->y;
		//printf("Coordinate=(%d,%d)\n",x[i],y[i]);
	}
	
	for(int i=0;i<length;i++)
	{
		sum_x=0,sum_y=0,sum_xx=0,sum_yy=0,sum_xy=0,sum_xxx=0,sum_yyy=0,sum_xxy=0,sum_xyy=0;
		
		for(int j=-N;j<N;j++)
		{
			sum_x+=x[(i+j+length)%length];
			sum_y+=y[(i+j+length)%length];
			sum_xx+=x[(i+j+length)%length]*x[(i+j+length)%length];
			sum_yy+=y[(i+j+length)%length]*y[(i+j+length)%length];
			sum_xy+=x[(i+j+length)%length]*y[(i+j+length)%length];
			sum_xxx+=x[(i+j+length)%length]*x[(i+j+length)%length]*x[(i+j+length)%length];
			sum_yyy+=y[(i+j+length)%length]*y[(i+j+length)%length]*y[(i+j+length)%length];
			sum_xxy+=x[(i+j+length)%length]*x[(i+j+length)%length]*y[(i+j+length)%length];
			sum_xyy+=x[(i+j+length)%length]*y[(i+j+length)%length]*y[(i+j+length)%length];
		}
		C=(N*sum_xx-sum_x*sum_y);
		D=(N*sum_xy-sum_x*sum_y);
		E=N*sum_xxx+N*sum_xyy-(sum_xx+sum_yy)*sum_x;
		G=(N*sum_yy-sum_y*sum_y);
		H=N*sum_xxy+N*sum_yyy-(sum_xx+sum_yy)*sum_y;
		a=(H*D-E*G)/(C*G-D*D);
		b=(H*C-E*D)/(D*D-G*C);
		c=-(sum_xx+sum_yy+a*sum_x+b*sum_y)/N;
		R=sqrt(a*a+b*b-4*c)/2;
		
		temp->curvature=1/R;
		List_Replace(i,temp);
			
		///
		printf("\nNo. %d     R=%f",i,R);
	}/****/
	delete temp;
	temp=NULL;
	

}

void SuperContour::Contour_rotation(int N)
{
	if(N>=length)
	{
		printf("\nThe LENGTH of the contour is too small, won't calculate the curvature.\n");
		return;
	}
	if(Contourlist==NULL)
	{
		printf("Empty LIST.\n");
		return;
	}
	double sum_x,sum_y,sum_xx,sum_yy,sum_xy;
	double k1,k2,k;
	
	int *x=(int*)malloc(length*sizeof(int));
	int *y=(int*)malloc(length*sizeof(int));
	
	ContourNode* temp=new ContourNode;
	
	for(int i=0;i<length;i++)
	{
		List_Visit(i,temp);
		x[i]=temp->x;
		y[i]=temp->y;
	}
	
	for(int i=0;i<length;i++)
	{
		sum_x=0,sum_y=0,sum_xx=0,sum_yy=0,sum_xy=0;		
		for(int j=0;j<N;j++)
		{
			sum_x+=x[(i+j+length)%length];
			sum_y+=y[(i+j+length)%length];
			sum_xx+=x[(i+j+length)%length]*x[(i+j+length)%length];
			sum_yy+=y[(i+j+length)%length]*y[(i+j+length)%length];
			sum_xy+=x[(i+j+length)%length]*y[(i+j+length)%length];
		}
		k1=(N*sum_xy-sum_x*sum_y)/(N*sum_xx-sum_x*sum_x);
		sum_x=0,sum_y=0,sum_xx=0,sum_yy=0,sum_xy=0;
		for(int j=0;j<N;j++)
		{
			sum_x+=x[(i-j+length)%length];
			sum_y+=y[(i-j+length)%length];
			sum_xx+=x[(i-j+length)%length]*x[(i-j+length)%length];
			sum_yy+=y[(i-j+length)%length]*y[(i-j+length)%length];
			sum_xy+=x[(i-j+length)%length]*y[(i-j+length)%length];
		}
		k2=(N*sum_xy-sum_x*sum_y)/(N*sum_xx-sum_x*sum_x);
		
		k=(k1-k2)/(1+k1*k2);
		
		temp->rotation=k;
		List_Replace(i,temp);
			
		///
		printf("\nNo. %d     R=%f",i,k);
	}/****/
	delete temp;
	temp=NULL;
	

}

CvPoint SuperContour::Contour_Center()
{
	CvPoint center;
	double sum_x=0,sum_y=0;
	ContourNode* temp=new ContourNode;
	if(length<0)
	{
		printf("\nEmpty LENGTH.");
		return cvPoint(0,0);
	}
	sum_x=0;
	sum_y=0;
	for(int i=0;i<length;i++)
	{
		List_Visit(i,temp);
		sum_x+=pow(temp->x,1.0);
		sum_y+=pow(temp->y,1.0);
	}
	center=cvPoint(pow(sum_x/length,1),pow(sum_y/length,1));
	//center=cvPoint(sum_x/length,sum_y/length);
	
	printf("center=(%d,%d)",center.x,center.y);
	return center;
}

/***Definition for mathching functions***/
void SuperContour::Match_Sampling(int SpNumber)
{
	ContourNode* temp=new ContourNode;
	if(length<1)
	{
		printf("\nContour does not exist.");
		return;
	}
	else if(length<=SpNumber)
	{
		printf("Too LARGE sampling-number");
		return;
	}
	
	
	
	int lag=length/SpNumber;
	if(Match_Data==NULL)
	{
		Match_Data=new Data_Match;
		Match_Data->Sample_Point=(CvPoint*)malloc(SpNumber*sizeof(CvPoint));
	}
	Match_Data->Sample_Number=SpNumber;
	Match_Data->Shift_Point=Contour_Center();
	for(int i=0;i<SpNumber;i++)
	{
		List_Visit(i*lag,temp);
		Match_Data->Sample_Point[i].x=temp->x;
		Match_Data->Sample_Point[i].y=temp->y;
		printf("\n%d",(temp->x-Match_Data->Shift_Point.x)*(temp->x-Match_Data->Shift_Point.x)+(temp->y-Match_Data->Shift_Point.y)*(temp->y-Match_Data->Shift_Point.y));
		//printf("(%d,%d)\n",temp->x,temp->y);
	}
	
	
	
}

void SuperContour::Match_Angle(SuperContour Contour_Target,int SpNumber,double Rotate_angle)
{
	double sum_x=0,sum_y=0;
	double x1,x2,y1,y2;
	Match_Sampling(SpNumber);
	Contour_Target.Match_Sampling(SpNumber);
	for(int i=0;i<SpNumber;i++)
	{
		
		x1=cos(Rotate_angle/180*3.141592654)*(Contour_Target.Match_Data->Sample_Point[i].x-Contour_Target.Match_Data->Shift_Point.x)-sin(Rotate_angle/180*3.141592654)*(Contour_Target.Match_Data->Sample_Point[i].y-Contour_Target.Match_Data->Shift_Point.x);
		y1=sin(Rotate_angle/180*3.141592654)*(Contour_Target.Match_Data->Sample_Point[i].x-Contour_Target.Match_Data->Shift_Point.x)+cos(Rotate_angle/180*3.141592654)*(Contour_Target.Match_Data->Sample_Point[i].y-Contour_Target.Match_Data->Shift_Point.y);
		
		x2=Match_Data->Sample_Point[i].x-Match_Data->Shift_Point.x;
		y2=Match_Data->Sample_Point[i].y-Match_Data->Shift_Point.y;
		
		sum_x+=x1*x2;
		sum_y+=y1*y2;
		//sum_x+=(Match_Data->Sample_Point[i].x-Match_Data->Shift_Point.x)*(Contour_Target.Match_Data->Sample_Point[i].x-Contour_Target.Match_Data->Shift_Point.x);
		//sum_y+=(Match_Data->Sample_Point[i].y-Match_Data->Shift_Point.y)*(Contour_Target.Match_Data->Sample_Point[i].y-Contour_Target.Match_Data->Shift_Point.y);
	}
	printf("\n %f %f",Rotate_angle,sum_x+sum_y);
}/****/

void SuperContour::Match_Best(SuperContour Contour_Target,int SpNumber)
{
	double best=0,second=0;
	double sum_x=0,sum_y=0,sum=0,bestsum=0;
	double x1,x2,y1,y2;
	Match_Sampling(SpNumber);
	Contour_Target.Match_Sampling(SpNumber);
	for(double Rotate_angle=0;Rotate_angle<360;Rotate_angle+=5)
	{
		
		for(int i=0;i<SpNumber;i++)
		{
			
			x1=cos(Rotate_angle/180*3.141592654)*(Contour_Target.Match_Data->Sample_Point[i].x-Contour_Target.Match_Data->Shift_Point.x)-sin(Rotate_angle/180*3.141592654)*(Contour_Target.Match_Data->Sample_Point[i].y-Contour_Target.Match_Data->Shift_Point.x);
			y1=sin(Rotate_angle/180*3.141592654)*(Contour_Target.Match_Data->Sample_Point[i].x-Contour_Target.Match_Data->Shift_Point.x)+cos(Rotate_angle/180*3.141592654)*(Contour_Target.Match_Data->Sample_Point[i].y-Contour_Target.Match_Data->Shift_Point.y);
		
			x2=Match_Data->Sample_Point[i].x-Match_Data->Shift_Point.x;
			y2=Match_Data->Sample_Point[i].y-Match_Data->Shift_Point.y;
		
			sum_x+=x1*x2;
			sum_y+=y1*y2;
			//sum_x+=(Match_Data->Sample_Point[i].x-Match_Data->Shift_Point.x)*(Contour_Target.Match_Data->Sample_Point[i].x-Contour_Target.Match_Data->Shift_Point.x);
			//sum_y+=(Match_Data->Sample_Point[i].y-Match_Data->Shift_Point.y)*(Contour_Target.Match_Data->Sample_Point[i].y-Contour_Target.Match_Data->Shift_Point.y);
		}
		sum=sum_x+sum_y;
		if(sum>bestsum)
		{
			bestsum=sum;
			best=Rotate_angle;
		}
		
	}
	printf("\n %f",best);
}/****/

void SuperContour::Match_FitAxis()
{
	
	if(OriginContour==NULL)
	{
		printf("\nEmpty Contour");
		return;
	}
	cvFitLine(OriginContour, CV_DIST_L2,0,0.01,0.01,Match_Data->AxisParam);
	//printf("\nx0=%f y0=%f",Match_Data->AxisParam[0],Match_Data->AxisParam[1]);
	return;
}
void SuperContour::Contour_Big()
{
	int flag0=0,flag1=0;
	if(OriginContour==NULL)
	{
		printf("\nEmpty Contour");
		return;
	}
	center=Contour_Center();
	Match_FitAxis();
	
	
	Match_Data->Big[0]=cvPoint(0,0);
	Match_Data->Big[1]=cvPoint(0,0);
	
	float x0,y0,x1,y1,sumx0,sumy0,sumx1,sumy1;
	x0=Match_Data->AxisParam[0];
	y0=Match_Data->AxisParam[1];
	
	printf("\nx0=%f y0=%f",x0,y0);
	
	int count_p,count_n;
	double distance,dtemp;
	CvPoint* p=NULL;
	
	for(distance=0;distance<30000;distance+=1)
	{
		count_p=0;
		count_n=0;
		sumx0=0;
		sumy0=0;
		sumx1=0;
		sumy1=0;
		for(int i=0;i<length;i++)
		{
			p=(CvPoint*)cvGetSeqElem(OriginContour,i);
			x1=p->x-Match_Data->AxisParam[2];
			y1=p->y-Match_Data->AxisParam[3];
			dtemp=(x1*y0-x0*y1);
			//printf("\ndtemp=%d",dtemp);
			if(dtemp>distance)
			{
				//printf("\n+++");
				count_p++;
				sumx0+=p->x;
				sumy0+=p->y;
			}
			if(dtemp<-distance)
			{
				count_n++;
				sumx1+=p->x;
				sumy1+=p->y;
			}
		}
		if(count_p<=30&&count_p>0&&flag0==0)
		{
			flag0=1;
			Match_Data->Big[0].x=sumx0/count_p;
			Match_Data->Big[0].y=sumy0/count_p;
			printf("\n count_p=%d",count_p);
		}
		if(count_n<=30&&count_n>0&&flag1==0)
		{
			flag1=1;
			Match_Data->Big[1].x=sumx1/count_n;
			Match_Data->Big[1].y=sumy1/count_n;
			printf("\n count_n=%d",count_n);
		}
		//printf("\n(%d,%d)",Match_Data->Big[0].x,Match_Data->Big[0].y);
		//printf("\n(%d,%d)",Match_Data->Big[1].x,Match_Data->Big[1].y);
	}
}

void SuperContour::Transform(float cs, float sn, CvPoint shift, double times)
{
	
	if(OriginContour==NULL||length==0)
	{
		printf("\n Empty Contour.");
		return;
	}
	Contour_Center();
	Contour_Big();
	
	
	float x,y;
	CvPoint* p=NULL;
	for(int i=0;i<length;i++)
	{
		p=(CvPoint*)cvGetSeqElem(OriginContour,i);
		x=p->x;
		y=p->y;
		
		printf("\n%f,%f",x,y);
		p->x=times*(cs*(x-center.x)-sn*(y-center.y))+shift.x;
		p->y=times*(sn*(x-center.x)+cs*(y-center.y))+shift.y;
		
	}
}
void SuperContour::Normalize()
{
	printf("\n|||cs=%f,sn=%f|||\n",Match_Data->cs,Match_Data->sn);
	Transform(Match_Data->cs,-(Match_Data->sn),center, 1/Match_Data->times);
}

void SuperContour::Fish_NormParam(SuperContour *Template)
{
	if(Template->OriginContour==NULL)
	{
		printf("\nEmpty Template");
		return;
	}
	Template->Contour_Big();
	
	float x_target=0,y_target=0,x_template=0,y_template=0;
	float cs=1,sn=0;
	
	
	float x0,y0,x1,y1;
	
	CvPoint cent_target,cent_template;
	
	cent_target=center;
	cent_template=Template->center;
	
	x_target=(Match_Data->Big[0].x+Match_Data->Big[1].x)/2-cent_target.x;
	y_target=(Match_Data->Big[0].y+Match_Data->Big[1].y)/2-cent_target.y;
	
	x_template=(Template->Match_Data->Big[0].x+Template->Match_Data->Big[1].x)/2-cent_template.x;
	y_template=(Template->Match_Data->Big[0].y+Template->Match_Data->Big[1].y)/2-cent_template.y;
	
	x1=x_target/sqrt(x_target*x_target+y_target*y_target);		///Normalized Direction Vector
	y1=y_target/sqrt(x_target*x_target+y_target*y_target);
	
	x0=x_template/sqrt(x_template*x_template+y_template*y_template);
	y0=y_template/sqrt(x_template*x_template+y_template*y_template);
	
	cs=x0*x1+y0*y1;
	sn=-x1*y0+x0*y1;
	Match_Data->cs=cs;
	Match_Data->sn=sn;
	Match_Data->times=sqrt(double((Match_Data->Big[0].x-Match_Data->Big[1].x)*(Match_Data->Big[0].x-Match_Data->Big[1].x)+(Match_Data->Big[0].y-Match_Data->Big[1].y)*(Match_Data->Big[0].y-Match_Data->Big[1].y))/((Template->Match_Data->Big[0].x-Template->Match_Data->Big[1].x)*(Template->Match_Data->Big[0].x-Template->Match_Data->Big[1].x)+(Template->Match_Data->Big[0].y-Template->Match_Data->Big[1].y)*(Template->Match_Data->Big[0].y-Template->Match_Data->Big[1].y)));
	printf("\n times=%f \n",Match_Data->times);
	printf("\n!!cs=%f sn=%f!!\n",cs,sn);
}

void SuperContour::Contour_Distance()
{
	if(length==0||OriginContour==NULL)
	{
		printf("\nEmpty Contour");
		return;
	}
	Contour_Center();
	ContourNode* output=NULL;
	for(int i=0;i<length;i++)
	{
		output=List_Use(i);
		if(output!=NULL)
		{
			output->distance=(output->x-center.x)*(output->x-center.x)+(output->y-center.y)*(output->y-center.y);
			if(output->distance<Match_Data->min_dist)
			{
				Match_Data->min_dist=output->distance;
			}
		}
		//printf("\nNo.%d,     distance=%d",i,output->distance);
	}
	
}
void SuperContour::Match_Fin()
{
	if(length==0)
	{
		printf("Empty Contour");
		return;
	}
	Contour_Center();
	Contour_Distance();
	Match_FitAxis();
	//Check step 1: maximum or minimum
	ContourNode *p=NULL, *np=NULL, *bp=NULL, *nnp=NULL, *bbp=NULL;
	
	for(int i=0;i<length;i+=3)
	{
		p=List_Use((i+length)%length);
		np=List_Use((i+10+length)%length);
		nnp=List_Use((i+3+length)%length);
		
		bp=List_Use((i-3+length)%length);
		bbp=List_Use((i-10+length)%length);
		printf(" %d   ",p->number);
		//if(p->distance>(np->distance+nnp->distance)/2&&p->distance>(bp->distance+bbp->distance)/2)
		if(p->distance>np->distance&&p->distance>bp->distance&&p->distance>nnp->distance&&p->distance>bbp->distance)			
		{
			p->peak=1;
			printf("1 \n");
		}
		else if(p->distance<(np->distance+nnp->distance)/2&&p->distance<(bp->distance+bbp->distance)/2)
		{
			p->peak=-1;
			printf("-1 \n");
		}
		else
		{
			p->peak=0;
			printf("0 \n");
		}
	}
	
	ContourNode *ptest[3]={NULL};
	ContourNode candidate[10];
	CvPoint vect=cvPoint(0,0);
	CvPoint Axvect=cvPoint(20.0*(Match_Data->AxisParam[0]),20.0*(Match_Data->AxisParam[1]));
	printf("Direction Axi(%d,%d)",Axvect.x,Axvect.y);
	
	int num_cand=0;
	num_cand=0;
	//Check step 2: min+max+min
	for(int i=0;i<length;i+=3)
	{
		p=List_Use((i+length)%length);
		if(p==NULL)
			printf("ErEr");
		if(p->peak==1||p->peak==-1)
		{
			ptest[0]=ptest[1];
			ptest[1]=ptest[2];
			ptest[2]=p;
			printf("\nHeihei\n");
		
			if(ptest[0]==NULL||ptest[1]==NULL||ptest[2]==NULL)
			{
				printf("\nHaha\n");
				continue;
			}
			else if(ptest[0]->peak==-1&&ptest[1]->peak==1&&ptest[2]->peak==-1)
			{
				printf("BiuBiu");
			
				if(abs(ptest[0]->number-ptest[2]->number)<length/5&&ptest[1]->distance-Match_Data->min_dist<10000)
				{
					
					//printf("        location (%d,%d)       ",ptest[1]->x,ptest[1]->y);
					num_cand++;
					//candidate=(ContourNode*)realloc(candidate,num_cand);
					if(candidate==NULL)
						break;
					printf("    num_cand=%d        ",num_cand);
					*(candidate+num_cand-1)=**(ptest+1);
					//*(candidate+num_cand-1);
					//**(ptest+1);
					printf("  fin NO= %d",ptest[1]->number);
				}
			}
		}
	}
	//Check step 3: Left or right and other conditions.
	
	for(int i=0;i<num_cand;i++)
	{
		vect=cvPoint((candidate+i)->x-center.x,(candidate+i)->y-center.y);
		if(vect.x*Axvect.y-vect.y*Axvect.x>0)
		{
			Match_Data->Fin_Left=cvPoint((candidate+i)->x,(candidate+i)->y);
		}
		if(vect.x*Axvect.y-vect.y*Axvect.x<0)
		{
			Match_Data->Fin_Right=cvPoint((candidate+i)->x,(candidate+i)->y);
		}
	}
}

/**
void SuperContour::Match_Fin()
{
	if(length==0)
	{
		printf("Empty Contour");
		return;
	}
	
	ContourNode *p=NULL, *np=NULL, *bp=NULL;
	int distance[20];
	
	double n,sx,sx2,sx3,sx4,sy,sxy,sx2y;
	double a,b,c;
	
	Contour_Distance();
	for(int i=0;i<length;i+=1)
	{
		n=20;
		sx=sx2=sx3=sx4=0;
		sy=sxy=sx2y=0;		
		a=b=c=0;
		
		for(int j=0;j<20;j++)
		{
			p=List_Use((i+j+length)%length);
			distance[j]=p->distance;
			//printf("\n%d",distance[j]);
		}
		for(int j=-10;j<10;j++)
		{
			sx+=j;
			sx2+=j*j;
			sx3+=j*j*j;
			sx4+=j*j*j*j;
			
			sy+=distance[(i+j+length)%length];
			sxy+=j*distance[(i+j+length)%length];
			sx2y+=j*j*distance[(i+j+length)%length];
		}
		a=(sx2*sx2*sx2y-sx*sx2y*sx3+sx*sx4*sxy+sx3*sx3*sy-sx2*(sx3*sxy+sx4*sy))/(sx2*sx2*sx2+n*sx3*sx3+sx*sx*sx4-sx2*(2*sx*sx3+n*sx4));
		b=(-sx*sx2*sx2y+n*sx2y*sx3+sx2*sx2*sxy-n*sx4*sxy-sx2*sx3*sy+sx*sx4*sy)/(sx2*sx2*sx2+n*sx3*sx3+sx*sx*sx4-sx2*(2*sx*sx3+n*sx4));
		c=(sx*sx*sx2y-n*sx2*sx2y+n*sx3*sxy+sx2*sx2*sy-sx*(sx2*sxy+sx3*sy))/(sx2*sx2*sx2+n*sx3*sx3+sx*sx*sx4-sx2*(2*sx*sx3+n*sx4));
		if(-b/2/a>0&&-b/2/a<20)
			printf("\n 1");
		else
			printf("\n 0");
		
	}
}**/