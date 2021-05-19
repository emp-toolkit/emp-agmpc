#include <emp-tool/emp-tool.h>
#include "emp-agmpc/emp-agmpc.h"
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
	FpreMP<nP> mp(ios, &pool, party);

	int num_ands = 1<<15;
	block * mac[nP+1];
	block * key[nP+1];
	bool * value;

	for(int i = 1; i <= nP; ++i) {
		key[i] = new block[num_ands*3];
		mac[i] = new block[num_ands*3];
	}
	value = new bool[num_ands*3];
	auto t1 = clock_start();
	mp.compute(mac, key, value, num_ands);
	cout <<"Gates: "<<num_ands<<" time: "<<time_from(t1)<<endl;
	return 0;
}
