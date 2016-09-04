#include <windows.h>

class UsClock{
	private:
		LARGE_INTEGER nFreq, nBeginTime, nEndTime;
	public:
		UsClock(){
			QueryPerformanceFrequency(&nFreq);
		}
		void setStartTime(){
			QueryPerformanceCounter(&nBeginTime);
		}
		void setEndTime(){
			QueryPerformanceCounter(&nEndTime);
		}
		double getTime(){
			return (double)(nEndTime.QuadPart-nBeginTime.QuadPart)/(double)nFreq.QuadPart;
		}
};
