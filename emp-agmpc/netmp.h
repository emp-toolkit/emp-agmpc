#ifndef NETIOMP_H__
#define NETIOMP_H__
#include <emp-tool/emp-tool.h>
#include "cmpc_config.h"
using namespace emp;

template<int nP>
class NetIOMP { public:
	NetIO*ios[nP+1];
	NetIO*ios2[nP+1];
	int party;
	bool sent[nP+1];
	NetIOMP(int party, int port) {
		this->party = party;
		memset(sent, false, nP+1);
		for(int i = 1; i <= nP; ++i)for(int j = 1; j <= nP; ++j)if(i < j){
			if(i == party) {
#ifdef LOCALHOST
				ios[j] = new NetIO(IP[j], port+2*(i*nP+j), true);
#else
				ios[j] = new NetIO(IP[j], port+2*(i), true);
#endif
				ios[j]->set_nodelay();	

#ifdef LOCALHOST
				ios2[j] = new NetIO(nullptr, port+2*(i*nP+j)+1, true);
#else
				ios2[j] = new NetIO(nullptr, port+2*(j)+1, true);
#endif
				ios2[j]->set_nodelay();	
			} else if(j == party) {
#ifdef LOCALHOST
				ios[i] = new NetIO(nullptr, port+2*(i*nP+j), true);
#else
				ios[i] = new NetIO(nullptr, port+2*(i), true);
#endif
				ios[i]->set_nodelay();	

#ifdef LOCALHOST
				ios2[i] = new NetIO(IP[i], port+2*(i*nP+j)+1, true);
#else
				ios2[i] = new NetIO(IP[i], port+2*(j)+1, true);
#endif
				ios2[i]->set_nodelay();	
			}
		}
	}
	int64_t count() {
		int64_t res = 0;
		for(int i = 1; i <= nP; ++i) if(i != party){
			res += ios[i]->counter;
			res += ios2[i]->counter;
		}
		return res;
	}

	~NetIOMP() {
		for(int i = 1; i <= nP; ++i)
			if(i != party) {
				delete ios[i];
				delete ios2[i];
			}
	}
	void send_data(int dst, const void * data, size_t len) {
		if(dst != 0 and dst!= party) {
			if(party < dst)
				ios[dst]->send_data(data, len);
			else
				ios2[dst]->send_data(data, len);
			sent[dst] = true;
		}
#ifdef __MORE_FLUSH
		flush(dst);
#endif
	}
	void recv_data(int src, void * data, size_t len) {
		if(src != 0 and src!= party) {
			if(sent[src])flush(src);
			if(src < party)
				ios[src]->recv_data(data, len);
			else
				ios2[src]->recv_data(data, len);
		}
	}
	NetIO*& get(size_t idx, bool b = false){
		if (b)
			return ios[idx];
		else return ios2[idx];
	}
	void flush(int idx = 0) {
		if(idx == 0) {
			for(int i = 1; i <= nP; ++i)
				if(i != party) {
					ios[i]->flush();
					ios2[i]->flush();
				}
		} else {
			if(party < idx)
				ios[idx]->flush();
			else
				ios2[idx]->flush();
		}
	}
	void sync() {
		for(int i = 1; i <= nP; ++i) for(int j = 1; j <= nP; ++j) if(i < j) {
			if(i == party) {
				ios[j]->sync();
				ios2[j]->sync();
			} else if(j == party) {
				ios[i]->sync();
				ios2[i]->sync();
			}
		}
	}
};
#endif //NETIOMP_H__
