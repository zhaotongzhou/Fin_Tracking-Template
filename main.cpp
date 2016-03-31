#define NUMBER 50


#include <cv.h>
#include <highgui.h>
#include <time.h>
#include <stdio.h>
#include "SuperContour.cpp"

void LongestContour(CvSeq* contours, CvSeq** ContourOfInterest);

class station
{
	CvPoint points[NUMBER];
	public:
		bool status=0; //if pushed
		CvPoint* p=NULL;//A pointer for reading data
		
		void init();
		void update(CvPoint point);
		void del(int index);
		void read(int index);
		
		int distribute();
		void eatcontour(CvSeq* contour,int threshold,CvPoint* fin);
		void returnfin(CvPoint* fin);
};

void station::init()
{
	for(int i=0;i<NUMBER;i++)
		points[i]=cvPoint(0,0);
	status=0;
}

void station::update(CvPoint point)
{
	if(status==0);
		status=1;
	
	for(int i=0;i<NUMBER-1;i++)
	{
		points[i]=points[i+1];
	}
	points[NUMBER-1]=cvPoint(point.x,point.y);
}

int station::distribute()
{
	int sumx=0,sumy=0;
	int sumxx=0,sumyy=0;
	int sum=0;
	
	if(status==0)
		return -10;
	for(int i=0;i<NUMBER;i++)
	{
		sumx+=points[i].x;
		sumy+=points[i].y;
		
		sumxx+=points[i].x*points[i].x;
		sumyy+=points[i].y*points[i].y;
	}
	sumx=sumx/NUMBER;
	sumy=sumy/NUMBER;
	
	sumxx/=NUMBER;
	sumyy/=NUMBER;
	
	sum=(sumxx+sumyy)-(sumx*sumx+sumy*sumy);
	return sum;
}

void station::eatcontour(CvSeq* contour,int threshold,CvPoint* fin)
{
	if(fin==NULL)
		fin=(CvPoint*)malloc(1*sizeof(CvPoint));
	int col=0;
	station::init();
	for(int i=0;i<contour->total;i++)
	{
		p=(CvPoint*)cvGetSeqElem(contour,i);
		
		update(*p);
		if(i%NUMBER==NUMBER-1)
		{
			col=distribute();
			
			
			//fin[i/NUMBER]=cvPoint(0,0);
			if(col<threshold&&fin!=NULL)
			{
				printf("\n pp %d\n",fin==NULL );
				fin=(CvPoint*)realloc(fin,((i/NUMBER)+1)*sizeof(CvPoint));
				//fin=(CvPoint*)realloc(fin,sizeof(CvPoint));
				
				printf("i/NUMBER=%d ",i/NUMBER);
				returnfin(fin+int(i/NUMBER));
				//returnfin(fin);
			}
			//printf("\nclo=%d\n",col);
		}
	}
}

void station::returnfin(CvPoint* fin)
{
	if(fin==NULL)
	{
		printf("Error!");
		return;
	}
	int sumx=0,sumy=0;
	for(int i=0;i<NUMBER;i++)
	{
		sumx+=points[i].x;
		sumy+=points[i].y;
	}
	fin->x=sumx/NUMBER;
	fin->y=sumy/NUMBER;
}


int main(int argc, char** argv)
{
	cvNamedWindow("tt",CV_WINDOW_AUTOSIZE);
	CvMemStorage* storage=cvCreateMemStorage(0);
	CvSeq* contour=NULL;
	CvSeq* fcontour=NULL;
	CvCapture* capture=NULL;
	
	capture  =cvCreateFileCapture(argv[1]);
	
	IplImage* frame;
	
	station checkstation;///This is for checking the distribution.
	CvPoint* fin=NULL;
	
	SuperContour fish;
	
	
	while(1)
	{
		fish.Super_Init();
		
		frame=cvQueryFrame(capture);
		if(!frame)
			break;
		IplImage* gray_img = cvCreateImage(cvGetSize(frame), IPL_DEPTH_8U, 1);
		cvCvtColor(frame, gray_img, CV_BGR2GRAY);
		//cvCopy(frame,gray_img);
		cvThreshold(gray_img, gray_img, 150, 150, CV_THRESH_BINARY);
		cvFindContours(gray_img,storage,&contour,sizeof(CvContour),CV_RETR_EXTERNAL);
		LongestContour(contour,&fcontour);
		
		cvDrawContours(gray_img,fcontour,cvScalarAll(255),cvScalarAll(255),0);
		//checkstation.eatcontour(contour,700,fin);
		//cvCircle(gray_img,*fin,4,cvScalarAll(255),3,8,0);
		
		fish.List_Create(fcontour);
		fish.Contour_Center();
		fish.Contour_Distance();
		fish.Match_FitAxis();
		
		fish.Match_Fin();
		cvCircle(gray_img,fish.Match_Data->Fin_Left,10,cvScalarAll(200));
		cvCircle(gray_img,fish.Match_Data->Fin_Right,10,cvScalarAll(200));
		cvCircle(gray_img,fish.center,5,cvScalarAll(200));
		
		
		fish.List_Clear();
		cvShowImage("tt",gray_img);
		cvZero(gray_img);
		char c=cvWaitKey(33);
		if(c==27) break;
	}
	cvReleaseCapture(&capture);
	cvDestroyWindow("tt");
	/****/
	
}

void LongestContour(CvSeq* contours, CvSeq** ContourOfInterest){
	CvSeq* biggestContour;
	//printf("---Finding Longest Contour---\n");
	int biggest=0;
		for (contours; contours!=NULL; contours=contours->h_next){
		//printf("%d elements\n",contours->total);
		if (contours->total > biggest){
			biggest=contours->total;
			biggestContour=contours;
			//printf("Currently the biggest!\n");
		}
	}

	*ContourOfInterest=cvCloneSeq(biggestContour);
}