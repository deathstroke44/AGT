#include "Preprocess.h"
#include "basis.h"
#include <fstream>
#include <assert.h>
#include <random>
#include <iostream>
#include <fstream>
#include <map>
#include <ctime>
#include <sstream>
#include <numeric>
<<<<<<< HEAD:cppCode/LSH-APG/src/Preprocess.cpp
#include <algorithm>
#include <omp.h>
=======
#include<algorithm>
using namespace std; 


>>>>>>> efb3637 (make changes for our experiments paradigms):dbLSH/src/Preprocess.cpp
#define E 2.718281746
#define PI 3.1415926

#define min(a,b)            (((a) < (b)) ? (a) : (b))
#define max(a,b)            (((a) > (b)) ? (a) : (b))

Preprocess::Preprocess(int _k, const std::string& path, const std::string& ben_file_)
{
	lsh::timer timer;
	k = _k;
	std::cout << "LOADING DATA..." << std::endl;
	timer.restart();
	load_data(path);
	std::cout << "LOADING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

	data_file = path;
	ben_file = ben_file_;
	if (data.N > 500) {
		ben_create();
	}
}

Preprocess::Preprocess(const std::string& path, const std::string& ben_file_, float beta_)
{

	hasT = true;
	beta = beta_;
	lsh::timer timer;
	std::cout << "LOADING DATA..." << std::endl;
	timer.restart();
	load_data(path);
	std::cout << "LOADING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

	data_file = path;
	ben_file = ben_file_;
	if (data.N > 1000) {
		ben_create();
	}
}


vector<vector<float>> readFVecsFromExternal(string filepath, int maxRow=-1) {
  FILE *infile = fopen(filepath.c_str(), "rb");
  vector<vector<float>> data= {};
  if (infile == NULL) {
    std::cout << "File not found" << std::endl;
    return data;
  }
  
  int rowCt = 0;
  int dimen;
  while (true) {
    if (fread(&dimen, sizeof(int), 1, infile) == 0) {
      break;
    }
    std::vector<float> v(dimen);
    if(fread(v.data(), sizeof(float), dimen, infile) == 0) {
      std::cout << "Error when reading" << std::endl;
    };
    
	data.push_back(v);

    rowCt++;
    
    if (maxRow != -1 && rowCt >= maxRow) {
      break;
    }
  }
  // std::cout<<"Row count test: "<<rowCt<<std::endl;

  if (fclose(infile)) {
    std::cout << "Could not close data file" << std::endl;
  }
  return data;
}

void Preprocess::load_data(const std::string& path)
{
<<<<<<< HEAD:cppCode/LSH-APG/src/Preprocess.cpp
	std::string file = path + "_new";
	std::ifstream in(file.c_str(), std::ios::binary);
	while (!in) {
		printf("Fail to find data file!\n");
		exit(0);
	}

	unsigned int header[3] = {};
	assert(sizeof header == 3 * 4);
	in.read((char*)header, sizeof(header));
	assert(header[0] == sizeof(float));
	data.N = header[1];
	data.dim = header[2];

	data.val = new float* [data.N];
	for (int i = 0; i < data.N; ++i) {
		data.val[i] = new float[data.dim];
		//in.seekg(sizeof(float), std::ios::cur);
		in.read((char*)data.val[i], sizeof(float) * header[2]);
	}

<<<<<<< HEAD:cppCode/LSH-APG/src/Preprocess.cpp
	//data.val = new float* [data.N];
	//float* dataBase = new float[data.N * data.dim];
	//in.read((char*)dataBase, sizeof(float) * (size_t)data.N * data.dim);
	//for (int i = 0; i < data.N; ++i) {
	//	data.val[i] = dataBase + i * data.dim;
	//	//in.seekg(sizeof(float), std::ios::cur);
	//	//in.read((char*)data.val[i], sizeof(float) * header[2]);
	//}
	int MaxQueryNum = 200;
	data.query = data.val;
	data.val = &(data.query[MaxQueryNum]);
	data.N -= MaxQueryNum;
=======
=======
>>>>>>> a908c48 (Fixed some bugs):dbLSH/src/Preprocess.cpp
	vector<vector<float>> base = readFVecsFromExternal(path+"base.fvecs");
	data.N = base.size();
	data.dim = base[0].size();
	data.val = new float* [data.N];
	for (int i = 0; i < data.N; ++i) {
		data.val[i] = new float[data.dim];
		//in.seekg(sizeof(float), std::ios::cur);
		for (int j=0;j<base[i].size();j++) {
			data.val[i][j]=base[i][j];
		}
	}
	vector<vector<float>> query = readFVecsFromExternal(path+"query.fvecs");
	data.numQuery = query.size();
	data.query = new float* [data.numQuery];
	for (int i = 0; i < data.numQuery; ++i) {
		data.query[i] = new float[data.dim];
		for (int j=0;j<query[i].size();j++) {
			data.val[i][j]=query[i][j];
		}
	}
>>>>>>> efb3637 (make changes for our experiments paradigms):dbLSH/src/Preprocess.cpp

	std::cout << "Load from new file: " << path << "\n";
	std::cout << "N=    " << data.N << "\n";
	std::cout << "dim=  " << data.dim << "\n\n";
	std::cout << "Qnum=  " << data.numQuery << "\n\n";
	std::cout << "size of base=  " << sizeof(&(data.val)) << "\n\n";
	std::cout << "size of query=  " << sizeof(&(data.query)) << "\n\n";
}

struct Tuple
{
	unsigned id;
	float dist;
};

bool comp(const Tuple& a, const Tuple& b)
{
	return a.dist < b.dist;
}

void Preprocess::ben_make()
{
<<<<<<< HEAD:cppCode/LSH-APG/src/Preprocess.cpp
	int MaxQueryNum = min(200, (int)data.N - 201);
	benchmark.N = MaxQueryNum, benchmark.num = 100;

	benchmark.N += 200;

=======
	int MaxQueryNum = data.numQuery;
	benchmark.N = MaxQueryNum, benchmark.num = k;
	cout<<"Benchmark k updated to: "<<k<<endl;
>>>>>>> efb3637 (make changes for our experiments paradigms):dbLSH/src/Preprocess.cpp
	benchmark.indice = new int* [benchmark.N];
	benchmark.dist = new float* [benchmark.N];
	for (unsigned j = 0; j < benchmark.N; j++) {
		benchmark.indice[j] = new int[benchmark.num];
		benchmark.dist[j] = new float[benchmark.num];
	}
	
	lsh::progress_display pd(benchmark.N);

#pragma omp parallel for num_threads(100)
	for (int j = 0; j < benchmark.N; j++)
	{
		std::vector<Tuple> dists(data.N);

		for (unsigned i = 0; i < data.N; i++)
		{
			dists[i].id = i;
			dists[i].dist = cal_distSqrt(data.val[i], data.query[j], data.dim);
		}

		sort(dists.begin(), dists.end(), comp);
		for (unsigned i = 0; i < benchmark.num; i++)
		{
			benchmark.indice[j][i] = (int)dists[i].id;
			benchmark.dist[j][i] = dists[i].dist;
		}
		++pd;
	}
}

void Preprocess::ben_save()
{
	std::ofstream out(ben_file.c_str(), std::ios::binary);
	out.write((char*)&benchmark.N, sizeof(unsigned));
	out.write((char*)&benchmark.num, sizeof(unsigned));

	for (unsigned j = 0; j < benchmark.N; j++) {
		out.write((char*)&benchmark.indice[j][0], sizeof(int) * benchmark.num);
	}

	for (unsigned j = 0; j < benchmark.N; j++) {
		out.write((char*)&benchmark.dist[j][0], sizeof(float) * benchmark.num);
	}

	out.close();
}

//Be cautious to use!!!!!!!!!!!!!!!!
//Please delete the sentence that uses this function after using it!!!!
void Preprocess::ben_correct()
{
	std::ofstream out(ben_file.c_str(), std::ios::binary);
	out.write((char*)&benchmark.N, sizeof(unsigned));
	out.write((char*)&benchmark.num, sizeof(unsigned));

	for (unsigned j = 0; j < benchmark.N; j++) {
		out.write((char*)&benchmark.indice[j][0], sizeof(int) * benchmark.num);
	}

	for (unsigned j = 0; j < benchmark.N; j++) {
		for (int l = 0; l < benchmark.num; ++l) {
			benchmark.dist[j][l] = sqrt(benchmark.dist[j][l]);
		}
		out.write((char*)&benchmark.dist[j][0], sizeof(float) * benchmark.num);
	}

	out.close();
	exit(0);
}

//Be cautious to use!!!!!!!!!!!!!!!!
//Please delete the sentence that uses this function after using it!!!!
void Preprocess::ben_correct_inverse()
{
	std::ofstream out(ben_file.c_str(), std::ios::binary);
	out.write((char*)&benchmark.N, sizeof(unsigned));
	out.write((char*)&benchmark.num, sizeof(unsigned));

	for (unsigned j = 0; j < benchmark.N; j++) {
		out.write((char*)&benchmark.indice[j][0], sizeof(int) * benchmark.num);
	}

	for (unsigned j = 0; j < benchmark.N; j++) {
		for (int l = 0; l < benchmark.num; ++l) {
			benchmark.dist[j][l] = (benchmark.dist[j][l]) * (benchmark.dist[j][l]);
		}
		out.write((char*)&benchmark.dist[j][0], sizeof(float) * benchmark.num);
	}

	out.close();
	exit(0);
}

void Preprocess::ben_load()
{
	std::ifstream in(ben_file.c_str(), std::ios::binary);
	in.read((char*)&benchmark.N, sizeof(unsigned));
	in.read((char*)&benchmark.num, sizeof(unsigned));

	benchmark.indice = new int* [benchmark.N];
	benchmark.dist = new float* [benchmark.N];
	for (unsigned j = 0; j < benchmark.N; j++) {
		benchmark.indice[j] = new int[benchmark.num];
		in.read((char*)&benchmark.indice[j][0], sizeof(int) * benchmark.num);
	}

	for (unsigned j = 0; j < benchmark.N; j++) {
		benchmark.dist[j] = new float[benchmark.num];
		in.read((char*)&benchmark.dist[j][0], sizeof(float) * benchmark.num);
	}
	in.close();
}

void Preprocess::ben_create()
{
	//unsigned a_test = data.N + 1;
	lsh::timer timer;
	std::ifstream in(ben_file.c_str(), std::ios::binary);
	//in.read((char*)&a_test, sizeof(unsigned));
	bool f = in.good();
	in.close();
<<<<<<< HEAD:cppCode/LSH-APG/src/Preprocess.cpp
	if (f){
=======
	if (a_test > 0 && a_test < data.N)//�ж��Ƿ��ܸ�дa_test
	{
>>>>>>> efb3637 (make changes for our experiments paradigms):dbLSH/src/Preprocess.cpp
		std::cout << "LOADING BENMARK..." << std::endl;
		timer.restart();
		ben_load();
		std::cout << "LOADING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

	}
	else{
		std::cout << "MAKING BENMARK..." << std::endl;
		timer.restart();
		ben_make();
		std::cout << "MAKING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

		std::cout << "SAVING BENMARK..." << std::endl;
		timer.restart();
		ben_save();
		std::cout << "SAVING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

	}
}

void Preprocess::showDataset()
{
	std::vector<float> lens(data.N);
	float sqrSum = 0, vSum = 0;
	float vmin = FLT_MAX, vmax = -FLT_MAX;
	float lmin = vmin, lmax = vmax;
	for (int i = 0; i < data.N; ++i) {
		float distSqr = cal_lengthSquare(data.val[i], data.dim);
		lens[i] = sqrt(distSqr);
		sqrSum += distSqr;
		lmin = min(lmin, lens[i]);
		lmax = max(lmax, lens[i]);
		for (int j = 0; j < data.dim; ++j) {
			vSum += data.val[i][j];
			vmin = min(vmin, data.val[i][j]);
			vmax = max(vmax, data.val[i][j]);
		}
	}

	float vavg = vSum / (data.N * data.dim);
	float vstd = (sqrSum) / (data.N * data.dim) - pow(vavg, 2);
	float lSum = std::accumulate(lens.begin(), lens.end(), 0.0f);
	float lavg = lSum / data.N;
	float lstd = (sqrSum) / (data.N) - -pow(lavg, 2);

	std::cout << "min value:  " << vmin << "\n";
	std::cout << "max value:  " << vmax << "\n";
	std::cout << "std value:  " << vstd << "\n";
	std::cout << "avg value:  " << vavg << "\n";
	//std::cout << "nstd value: " << vavg << "\n";
	std::cout << "min dist:   " << lmin << "\n";
	std::cout << "max dist:   " << lmax << "\n";
	std::cout << "std dist:   " << lstd << "\n";
	std::cout << "avg dist:   " << lavg << "\n";
	//std::cout << "nstd value: " << vavg << "\n";
	exit(0);
}

Preprocess::~Preprocess()
{
<<<<<<< HEAD:cppCode/LSH-APG/src/Preprocess.cpp
	int MaxQueryNum = min(200, (int)data.N - 201);
	clear_2d_array(data.query, data.N + MaxQueryNum);
=======
	int MaxQueryNum = data.numQuery;

	clear_2d_array(data.query, data.N + data.numQuery);
>>>>>>> efb3637 (make changes for our experiments paradigms):dbLSH/src/Preprocess.cpp
	//clear_2d_array(Dists, MaxQueryNum);
	clear_2d_array(benchmark.indice, benchmark.N);
	clear_2d_array(benchmark.dist, benchmark.N);
	delete[] SquareLen;
}


Parameter::Parameter(Preprocess& prep, unsigned L_, unsigned K_, float rmin_)
{
	N = prep.data.N;
	dim = prep.data.dim;
	L = L_;
	K = K_;
	MaxSize = 5;
	R_min = rmin_;
}

Parameter::~Parameter()
{
}