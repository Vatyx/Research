#include "mmap_file.h"

MMapFile mmap_file_read(string filename)
{
	MMapFile f;

	f.hfile = CreateFileA(filename.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (!f.hfile)
	{
		printf("ERROR: Failed to open file %s\n", filename.c_str());
		exit(1);
	}
	f.size = GetFileSize(f.hfile, 0); // TODO: make 64 bit
	f.hmmfile = CreateFileMapping (f.hfile, NULL, PAGE_READONLY, 0, 0, NULL);
	if (!f.hmmfile)
	{
		printf("ERROR: Failed to open file %s\n", filename.c_str());
		exit(1);
	}
	f.ptr = (unsigned char*) MapViewOfFile(f.hmmfile, FILE_MAP_READ, 0, 0, 0); 
	if (!f.ptr)
	{
		printf("ERROR: Failed to open file %s\n", filename.c_str());
		exit(1);
	}

	return f;
}

MMapFile mmap_file_write(string filename, __int64 size)
{
	MMapFile f;

	f.size = size;

	f.hfile = CreateFileA(filename.c_str(), GENERIC_READ | GENERIC_WRITE, 0, NULL, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (!f.hfile)
	{
		printf("ERROR: Failed to open file %s\n", filename.c_str());
		exit(1);
	}
	f.hmmfile = CreateFileMapping (f.hfile, NULL, PAGE_READWRITE, (__int32)(f.size >> 32), (__int32)(f.size & 0xffffffff), NULL);
	if (!f.hmmfile)
	{
		printf("ERROR: Failed to open file %s\n", filename.c_str());
		exit(1);
	}
	f.ptr = (unsigned char*) MapViewOfFile(f.hmmfile, FILE_MAP_ALL_ACCESS, 0, 0, 0); 
	if (!f.ptr)
	{
		printf("ERROR: Failed to open file %s\n", filename.c_str());
		exit(1);
	}

	return f;
}

void mmap_file_close(MMapFile &f)
{
	UnmapViewOfFile(f.ptr);
	CloseHandle(f.hmmfile);
	CloseHandle(f.hfile);
}

void load_file_to_string(string filename, string &contents)
{
	MMapFile f = mmap_file_read(filename);
	contents.resize((int)f.size);
	for (int i = 0; i < f.size; i++)
		contents[i] = f.ptr[i];
	mmap_file_close(f);
}
