#include <emp-tool/emp-tool.h>
#include "emp-agmpc/emp-agmpc.h"
using namespace std;
using namespace emp;

//const string circuit_file_location = macro_xstr(EMP_CIRCUIT_PATH);

// test from root directory of repo using
// ./run ./bin/pr10_test ...

const static int nP = 3;
int party, port;
void bench_once(NetIOMP<nP> * ios[2], ThreadPool * pool, string filename) {
	if(party == 1)cout <<"CIRCUIT:\t"<<filename<<endl;
	//string file = circuit_file_location+"/"+filename;
	BristolFormat cf(filename.c_str());

	auto start = clock_start();
	CMPC<nP>* mpc = new CMPC<nP>(ios, pool, party, &cf);
	ios[0]->flush();
	ios[1]->flush();
	double t2 = time_from(start);
//	ios[0]->sync();
//	ios[1]->sync();
	if(party == 1)cout <<"Setup:\t"<<party<<"\t"<< t2 <<"\n"<<flush;

	start = clock_start();
	mpc->function_independent();
	ios[0]->flush();
	ios[1]->flush();
	t2 = time_from(start);
	if(party == 1)cout <<"FUNC_IND:\t"<<party<<"\t"<<t2<<" \n"<<flush;

	start = clock_start();
	mpc->function_dependent();
	ios[0]->flush();
	ios[1]->flush();
	t2 = time_from(start);
	if(party == 1)cout <<"FUNC_DEP:\t"<<party<<"\t"<<t2<<" \n"<<flush;

	bool *in = new bool[cf.n1+cf.n2]; bool *out = new bool[cf.n3];
	memset(in, false, cf.n1+cf.n2);
	start = clock_start();
	mpc->online(in, out);
	ios[0]->flush();
	ios[1]->flush();
	t2 = time_from(start);
//	uint64_t band2 = io.count();
//	if(party == 1)cout <<"bandwidth\t"<<party<<"\t"<<band2<<endl;
	if(party == 1)cout <<"ONLINE:\t"<<party<<"\t"<<t2<<" \n"<<flush;
	delete mpc;
}
int main(int argc, char** argv) {
	parse_party_and_port(argv, &party, &port);
	if(party > nP)return 0;
	NetIOMP<nP> io(party, port);
#ifdef LOCALHOST
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
#else
	NetIOMP<nP> io2(party, port+2*(nP+1));
#endif
	NetIOMP<nP> *ios[2] = {&io, &io2};
	ThreadPool pool(2*(nP-1)+2);	

    char d[100];
    char* t = getcwd(d, 100);
    (void) t;
    string circuit_file_location = string(d) + "/test/";
	bench_once(ios, &pool, circuit_file_location+"pr10_circuit.txt");
	return 0;
}