#include <emp-tool>
#include "netmp.h"
#include "abitmp.h"
#include "fpremp.h"
using namespace std;

int main(int argc, char** argv) {
	int port, party;
	parse_party_and_port(argv, &party, &port);

	NetIOMP<3> io(party, port);
	for(int i = 1; i <= 3; ++i)for(int j = 1; j <= 3; ++j)if(i < j) {
		if(i == party) {
			int data = i*100+j;
			io.send_data(j, &data, 4);
		} else if (j == party) {
			int data = 0;
			io.recv_data(i, &data, 4);
			if(data!= i*100+j) {
				cout <<"WRONG\n"<<endl;
			}
		}
	}
	io.flush();
	ThreadPool pool(2*3*3);	
	block *MAC[4], *KEY[4];
	bool * data;
	PRG prg;
	int LEN = 1<<25;
	data = new bool[LEN*3];
	prg.random_bool(data, LEN*3);
	for(int i = 1; i <= 3 ; ++i) {
		MAC[i] = new block[LEN*3];
		KEY[i] = new block[LEN*3];
	}

	for(int LEN = 1<<10; LEN <= 1<<25; LEN*=2 ) {
		ABitMP<3> abitmp(&io, &pool, party);
		double t1 = timeStamp();
		abitmp.compute(MAC, KEY, data, LEN);
		if(party == 1)
			cout <<LEN<<"\t"<< (timeStamp()-t1)/(0.0+LEN)<<endl;
		auto ret = abitmp.check(MAC, KEY, data, LEN);
		ret.get();
	}
	delete[] data;
	for(int i = 1; i <= 3; ++i) {
		delete[] MAC[i];
		delete[] KEY[i];
	}

	///	FpreMP<3> fpre(&io, &pool, party);	
	//	fpre.compute(MAC, KEY, data, LEN);

	//	cout <<"DONE3 !"<<party<<" \n"<<flush;
	return 0;
}
