#pragma once

#include <string>
#include <windows.h>

using namespace std;

typedef unsigned char byte;

struct MMapFile
{
	HANDLE hfile;
	HANDLE hmmfile;
	byte *ptr;
	__int64 size;
};

MMapFile mmap_file_read(string filename);
MMapFile mmap_file_write(string filename, __int64 size);
void mmap_file_close(MMapFile &f);
void load_file_to_string(string filename, string &contents);
