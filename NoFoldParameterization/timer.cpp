#include "timer.h"
#include "stdio.h"

#ifdef _WIN32

#include <windows.h>
double get_time()
{
	static LARGE_INTEGER t_start;
	static bool first_time = true;
	if (first_time)
	{
		first_time = false;
		QueryPerformanceCounter(&t_start);
	}
	LARGE_INTEGER t, f;
	QueryPerformanceCounter(&t);
	QueryPerformanceFrequency(&f);
	double time = double(t.QuadPart - t_start.QuadPart)/double(f.QuadPart);
	return time;
}

#else

#include <sys/time.h>
#include <sys/resource.h>

double get_time()
{
	struct timeval t;
	struct timezone tzp;
	gettimeofday(&t, &tzp);
	return t.tv_sec + t.tv_usec*1e-6;
}

#endif
