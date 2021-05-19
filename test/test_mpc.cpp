#include <emp-tool/emp-tool.h>
#include "emp-agmpc/emp-agmpc.h"
using namespace std;
using namespace emp;

const string circuit_file_location = macro_xstr(EMP_CIRCUIT_PATH) + string("bristol_format/");
static char out3[] = "92b404e556588ced6c1acd4ebf053f6809f73a93";//bafbc2c87c33322603f38e06c3e0f79c1f1b1475";

int main(int argc, char** argv) {
	int port, party;
	parse_party_and_port(argv, &party, &port);

	const static int nP = 3;
	NetIOMP<nP> io(party, port);
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
	NetIOMP<nP> *ios[2] = {&io, &io2};
	ThreadPool pool(4);	
	string file = circuit_file_location+"/AES-non-expanded.txt";
	file = circuit_file_location+"/sha-1.txt";
	BristolFormat cf(file.c_str());

	CMPC<nP>* mpc = new CMPC<nP>(ios, &pool, party, &cf);
	cout <<"Setup:\t"<<party<<"\n";

	mpc->function_independent();
	cout <<"FUNC_IND:\t"<<party<<"\n";

	mpc->function_dependent();
	cout <<"FUNC_DEP:\t"<<party<<"\n";

	bool in[512]; bool out[160];
	memset(in, false, 512);
	mpc->online(in, out);
	uint64_t band2 = io.count();
	cout <<"bandwidth\t"<<party<<"\t"<<band2<<endl;
	cout <<"ONLINE:\t"<<party<<"\n";
	if(party == 1) {
		string res = "";
		for(int i = 0; i < cf.n3; ++i)
			res += (out[i]?"1":"0");
		cout << hex_to_binary(string(out3))<<endl;
		cout << res<<endl;
		cout << (res == hex_to_binary(string(out3))? "GOOD!":"BAD!")<<endl<<flush;
	}
	delete mpc;
	return 0;
}
