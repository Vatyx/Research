#include <stdlib.h>
#include "vect.h"
#include "array.h"
#include "obj_data.h"
#include "mmap_file.h"

#include "def_common.h"

int myatoi(char *s)
{
	int n = 0;
	while (*s != '\0') {
		n = n * 10 + *s - '0';
		s++;
	}

	return n;
}

double myatof(char *s)
{
	int mode = 0; // 0 = integer, 1 = decimal, 2 = exponent

	float n = 0;
	float t = 0.1f;
	int e = 0;
	int eq = 1;
	float q;

	if (*s == '-') {
		q = -1.0f;
		s++;
	}
	else
		q = 1.0f;

	for (;;)
	{
		if (mode == 0)
		{
			for (;;)
			{
				if (*s == '\0') 
					return q * n;
				else if (*s == '.')
				{
					mode = 1;
					s++;
					break;
				}
				else
				{
					n = n * 10 + *s - '0';
					s++;
				}
			}
		}
		else if (mode == 1)
		{
			for (;;)
			{
				if (*s == '\0') 
					return q * n;
				else if (*s == 'e')
				{
					mode = 2;
					s++;
					break;
				}
				else
				{
					n += (*s - '0') * t;
					t *= 0.1f;
					s++;
				}
			}
		}
		else if (mode == 2)
		{
			for (;;)
			{
				if (*s == '\0') 
					return q * n * pow(10.0, eq*e);
				else if (*s == '-')
				{
					eq = -1;
					s++;
				}
				else
				{
					e = e * 10 + *s - '0';
					s++;
				}
			}
		}
		else
			break;
	}

	return 0;
}

char* eat_line(char *ptr, char *end)
{
	for (;ptr < end; ptr++)
	{
		if (*ptr == '\n')
			break;
	}
	return ptr;
}

char* eat_white(char *ptr, char *end)
{
	for (;ptr < end; ptr++)
	{
		if (*ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r')
			break;
	}
	return ptr;
}

char* eat_space(char *ptr, char *end)
{
	for (;ptr < end; ptr++)
	{
		if (*ptr != ' ' && *ptr != '\t')
			break;
	}
	return ptr;
}

char* get_token(char *buf, char *ptr, char *end)
{
	ptr = eat_white(ptr, end);

	for (;ptr < end; ptr++, buf++)
	{
		*buf = *ptr;
		if (*ptr == ' ' || *ptr == '\t' || *ptr == '\n' || *ptr == '\r')
			break;
	}

	*buf = 0;
	return ptr;
}

char* get_obj_index(char *buf, char *ptr)
{
	// eat separater
	while (*ptr == '/')
		ptr++;

	for (;; ptr++, buf++)
	{
		*buf = *ptr;
		if (*ptr == '/' || *ptr == 0)
			break;
	}

	*buf = 0;
	return ptr;
}

#include <vector>
void load_obj(string filename, MeshData &mesh)
{
	vector<vect3d> pos;
	vector<vect2d> uv;
	vector<bool> found;
	// open the file
	MMapFile mmfile = mmap_file_read(filename);

	// read the file
	char *ptr = (char*)mmfile.ptr;
	char *end = (char*)(mmfile.ptr + mmfile.size);
	char buf[64];
//	int count = 0;
	for (;;)
	{
		ptr = get_token(buf, ptr, end);
		if (buf[0] == 'v' && buf[1] == 0)
		{
			vect3d v;

			for (int i = 0; i < 3; i++)
			{
				ptr = get_token(buf, ptr, end);
				v[i] = SCALE_MESH*(float)myatof(buf);
			}

			if(v[0] < -1000)
			{
				v[0] = 0;
			}

			if(v[1] < -1000)
			{
				v[1] = 0;
			}

			if(v[2] < -1000)
			{
				v[2] = 0;
			}

			if(v[0] < -1000)
			{
				printf("WTF\n");
			}
			pos.push_back(v);
			found.push_back(false);
		}
		else if (buf[0] == 'v' && buf[1] == 'n' && buf[2] == 0)
		{
			vect3d v;

			for (int i = 0; i < 3; i++)
			{
				ptr = get_token(buf, ptr, end);
				v[i] = (float)myatof(buf);
			}

			mesh.normals.push_back(v);
		}
		else if (buf[0] == 'v' && buf[1] == 't' && buf[2] == 0)
		{
			vect2d v;

			for (int i = 0; i < 2; i++)
			{
				ptr = get_token(buf, ptr, end);
				v[i] = SCALE_TEX_COORDS*(float)myatof(buf);
			}

			uv.push_back(v);
		}
		else if (buf[0] == 'f' && buf[1] == 0)
		{
			int idx_buf[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
			int idx_tex[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
			int idx_norm[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
			int idx_num = 0;
			do
			{
				ptr = get_token(buf, ptr, end);
				
				// parse the token for indices
				char *ptr_i = buf;
				char buf_i[64];
				
				ptr_i = get_obj_index(buf_i, ptr_i);
				idx_buf[idx_num] = myatoi(buf_i) - 1;

				ptr_i = get_obj_index(buf_i, ptr_i);
				idx_tex[idx_num] = myatoi(buf_i) - 1;
				
				ptr_i = get_obj_index(buf_i, ptr_i);
				idx_norm[idx_num] = myatoi(buf_i) - 1;

				idx_num++;
				assert(idx_num <= 8);
				ptr = eat_space(ptr, end);
			} while (*ptr != '\n' && *ptr != '\r' && ptr != end);

			for (int k = 2; k < idx_num; k++)
			{
				mesh.indices_pos.push_back(vect3i::make(idx_buf[0], idx_buf[k-1], idx_buf[k]));
				mesh.indices_tex.push_back(vect3i::make(idx_tex[0], idx_tex[k-1], idx_tex[k]));
				mesh.indices_norm.push_back(vect3i::make(idx_norm[0], idx_norm[k-1], idx_norm[k]));
			}
		}
		else
		{
			ptr = eat_line(ptr, end);
		}
		
		if (ptr >= end)
			break;
	}

	for (int i = 0; i < found.size(); i++)
	{
		vect3d t;
		t[0] = t[1] = t[2] = 0;

		vect2d tt;
		tt[0] = t[1] = 0;
		mesh.positions.push_back(t);
		mesh.tex_coords.push_back(tt);
	}
	for (int i = 0; i < mesh.indices_pos.s; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (!found[mesh.indices_tex[i][j]])
			{
				mesh.positions[mesh.indices_tex[i][j]] = pos[mesh.indices_pos[i][j]];
				mesh.tex_coords[mesh.indices_tex[i][j]] = uv[mesh.indices_tex[i][j]];
			}
		}
	}
	// close the file
	mmap_file_close(mmfile);
}

void load_obj2(string filename, MeshData &mesh)
{
	// open the file
	MMapFile mmfile = mmap_file_read(filename);

	// read the file
	char *ptr = (char*)mmfile.ptr;
	char *end = (char*)(mmfile.ptr + mmfile.size);
	char buf[64];
	for (;;)
	{
		ptr = get_token(buf, ptr, end);
		if (buf[0] == 'v' && buf[1] == 0)
		{
			vect3d v;

			for (int i = 0; i < 3; i++)
			{
				ptr = get_token(buf, ptr, end);
				v[i] = SCALE_MESH*(float)myatof(buf);
			}
			//v[0] *= -1;
			mesh.positions.push_back(v);

			vect2d vt;
			for (int i = 0; i < 2; i++)
			{
				ptr = get_token(buf, ptr, end);
				vt[i] = SCALE_TEX_COORDS*(float)myatof(buf);
			}

			//vt[0] *= -1;
			mesh.tex_coords.push_back(vt);
		}
		else if (buf[0] == 'v' && buf[1] == 'n' && buf[2] == 0)
		{
			printf("WHAT\n");

		}
		else if (buf[0] == 'v' && buf[1] == 't' && buf[2] == 0)
		{
			printf("WHAT\n");
		}
		else if (buf[0] == 'f' && buf[1] == 0)
		{
			int idx_buf[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
			int idx_tex[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
			int idx_norm[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
			int idx_num = 0;
			do
			{
				ptr = get_token(buf, ptr, end);
				
				// parse the token for indices
				char *ptr_i = buf;
				char buf_i[64];
				
				ptr_i = get_obj_index(buf_i, ptr_i);
				idx_buf[idx_num] = myatoi(buf_i) - 1;

				ptr_i = get_obj_index(buf_i, ptr_i);
				idx_tex[idx_num] = myatoi(buf_i) - 1;
				
				ptr_i = get_obj_index(buf_i, ptr_i);
				idx_norm[idx_num] = myatoi(buf_i) - 1;
//
				idx_num++;
				assert(idx_num <= 8);
			} while (*ptr != '\n' && *ptr != '\r');

			for (int k = 2; k < idx_num; k++)
			{
				mesh.indices_pos.push_back(vect3i::make(idx_buf[0], idx_buf[k-1], idx_buf[k]));
				mesh.indices_tex.push_back(vect3i::make(idx_buf[0], idx_buf[k-1], idx_buf[k]));
				mesh.indices_norm.push_back(vect3i::make(idx_norm[0], idx_norm[k-1], idx_norm[k]));
			}
		}
		else
		{
			ptr = eat_line(ptr, end);
		}
		
		if (ptr >= end)
			break;
	}

	// close the file
	mmap_file_close(mmfile);
}

void load_off(string filename, MeshData &mesh)
{
	// open the file
	MMapFile mmfile = mmap_file_read(filename);

	// read the file
	char *ptr = (char*)mmfile.ptr;
	char *end = (char*)(mmfile.ptr + mmfile.size);
	char buf[64];
	
	// read OFF header
	ptr = get_token(buf, ptr, end);

	ptr = get_token(buf, ptr, end);
	int num_verts = myatoi(buf);
	mesh.positions.resize(num_verts);

	ptr = get_token(buf, ptr, end);
	int num_faces = myatoi(buf);
	mesh.indices_pos.reserve(num_faces);
	
	ptr = get_token(buf, ptr, end);

	for (int ix = 0; ix < num_verts; ix++)
	{
		for (int i = 0; i < 3; i++)
		{
			ptr = get_token(buf, ptr, end);
			mesh.positions[ix][i] = (float)myatof(buf);
		}
	}

	for (int ix = 0; ix < num_faces; ix++)
	{
		int idx_buf[8];

		ptr = get_token(buf, ptr, end);
		int idx_num = myatoi(buf);
		for (int i = 0; i < idx_num; i++)
		{
			ptr = get_token(buf, ptr, end);
			idx_buf[i] = myatoi(buf);
		}

		for (int k = 2; k < idx_num; k++)
		{
			mesh.indices_pos.push_back(vect3i::make(idx_buf[0], idx_buf[k-1], idx_buf[k]));
		}
	}

	// close the file
	mmap_file_close(mmfile);
}
