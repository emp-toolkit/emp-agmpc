#include <emp-tool/emp-tool.h>
#include "emp-agmpc/emp-agmpc.h"
#include "test/test.h"
using namespace std;
using namespace emp;

const static int nP = 3;
int party, port;
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

	bench_once<nP>(party, ios, &pool, circuit_file_location+"sha-1.txt");
	return 0;
}
