
CXX = g++
CXXFLAGS= -c -v -Wall 
CXXFLAGS2= -c -v -Wall 


CVdir=C:/Program\ Files\ \(x86\)/OpenCV





#OpenCV Include directories (for header files)
openCVincludes = -I$(CVdir)/cxcore/include -I$(CVdir)/otherlibs/highgui -I$(CVdir)/cv/include

#openCVincludes =-IF:/opencv/opencv/build/include/
CVlibs= $(CVdir)/lib/cv.lib $(CVdir)/lib/highgui.lib $(CVdir)/lib/cxcore.lib
3rdpartyobjects= $(CVlibs)


all: main.exe

main.exe: main.o
	g++ -o main.exe main.o $(CVlibs)
		
main.o : main.cpp SuperContour.cpp
	$(CXX) $(CXXFLAGS2) -g main.cpp -o main.o $(openCVincludes)
#g++ -I/usr/include/opencv main.cpp -o main -lopencv_core -lopencv_imgproc -lopencv_highgui -lopencv_ml -lopencv_video -lopencv_features2d -lopencv_calib3d -lopencv_objdetect -lopencv_contrib -lopencv_legacy -lopencv_flann  



clean:

	 
	
	
