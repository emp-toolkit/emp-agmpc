#ifndef EMP_AGMPC_TEST_H__
#define EMP_AGMPC_TEST_H__
const string circuit_file_location = macro_xstr(EMP_CIRCUIT_PATH) + string("bristol_format/");


template<int nP>
int communication(NetIOMP<nP> * ios[2]) {
	return ios[0]->count() + ios[1]->count();
}

template<int nP>
void bench_once(int party, NetIOMP<nP> * ios[2], ThreadPool * pool, string filename) {
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
	if(party == 1) cout <<"FUNC_IND:\t"<<party<<"\t"<<t2<<" \n"<<flush;

	start = clock_start();
	mpc->function_dependent();
	ios[0]->flush();
	ios[1]->flush();
	t2 = time_from(start);
	if(party == 1) cout <<"FUNC_DEP:\t"<<party<<"\t"<<t2<<" \n"<<flush;
//	cout << "Offline Communication:\t "<<communication<nP>(ios)/1000.0/1000.0<<" MB"<<endl;

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
	cout << "Total Communication:\t "<<communication<nP>(ios)/1000.0/1000.0<<" MB"<<endl;
	delete mpc;
}

#endif// EMP_AGMPC_TEST_HHH