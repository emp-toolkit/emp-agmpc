#ifndef EMP_AGMPC_FLEXIBLE_INPUT_OUTPUT_H
#define EMP_AGMPC_FLEXIBLE_INPUT_OUTPUT_H

using namespace std;

template<int nP>
struct AuthBitShare
{
	bool bit_share;
	block key[nP + 1];
	block mac[nP + 1];
};

struct BitWithMac
{
	bool bit_share;
	block mac;
};

template<int nP>
class FlexIn
{
public:

	int len{};
	int party{};

	bool cmpc_associated = false;
	ThreadPool* pool;
	bool *value;
	block *key[nP + 1];
	block *mac[nP + 1];
	NetIOMP<nP> * io;
	block Delta;

	vector<int> party_assignment;
	// -2 represents an un-authenticated share (e.g., for random tape),
	// -1 represents an authenticated share,
	//  0 represents public input/output,

	vector<bool> plaintext_assignment; // if `party` provides the value for this bit, the plaintext value is here
	vector<AuthBitShare<nP>> authenticated_share_assignment; // if this bit is from authenticated shares, the authenticated share is stored here

	FlexIn(int len, int party) {
		this->len = len;
		this->party = party;
		this->pool = pool;

		AuthBitShare<nP> empty_abit;
		memset(&empty_abit, 0, sizeof(AuthBitShare<nP>));

		party_assignment.resize(len, 0);
		plaintext_assignment.resize(len, false);
		authenticated_share_assignment.resize(len, empty_abit);
	}

	~FlexIn() {
		party_assignment.clear();
		plaintext_assignment.clear();
		authenticated_share_assignment.clear();
	}

	void associate_cmpc(ThreadPool *associated_pool, bool *associated_value, block *associated_mac[nP + 1], block *associated_key[nP + 1],  NetIOMP<nP> *associated_io, block associated_Delta) {
		this->cmpc_associated = true;
		this->pool = associated_pool;
		this->value = associated_value;
		for(int j = 1; j <= nP; j++) {
			this->mac[j] = associated_mac[j];
			this->key[j] = associated_key[j];
		}
		this->io = associated_io;
		this->Delta = associated_Delta;
	}

	void assign_party(int pos, int which_party) {
		party_assignment[pos] = which_party;
	}

	void assign_plaintext_bit(int pos, bool cur_bit) {
		assert(party_assignment[pos] == party || party_assignment[pos] == -2  || party_assignment[pos] == 0);
		plaintext_assignment[pos] = cur_bit;
	}

	void assign_authenticated_bitshare(int pos, AuthBitShare<nP> *abit) {
		assert(party_assignment[pos] == -1);
		memcpy(&authenticated_share_assignment[pos], abit, sizeof(AuthBitShare<nP>));
	}

	void input(bool *masked_input_ret) {
		assert(cmpc_associated);

		/* assemble an array of the input masks, their macs, and their keys */
		/*
		 * Then,
		 * 		for a plaintext bit, the input mask, as well as its MAC, is sent to the input party, who uses the KEY for verification;
		 * 		for an un-authenticated bit, the input mask XOR with the input share is broadcast;
		 *		for an authenticated bit share, they are used to masked the previously data (and then checking its opening)
		 */
		vector<AuthBitShare<nP>> input_mask;
		for(int i = 0; i < len; i++) {
			AuthBitShare<nP> abit;
			abit.bit_share = value[i];
			for(int j = 1; j <= nP; j++) {
				if(j != party) {
					abit.key[j] = key[j][i];
					abit.mac[j] = mac[j][i];
				}
			}
			input_mask.emplace_back(abit);
		}

		/*
		 * first of all, handle the case party_assignment[] > 0
		 */

		/* prepare the bit shares to open for the corresponding party */
		vector<vector<BitWithMac>> open_bit_shares_for_plaintext_input_send;
		open_bit_shares_for_plaintext_input_send.resize(nP + 1);

		for(int j = 1; j <= nP; j++) {
			if(j != party) {
				for(int i = 0; i < len; i++) {
					BitWithMac mbit{};
					if(party_assignment[i] == j) {
						mbit.bit_share = input_mask[i].bit_share;
						mbit.mac = input_mask[i].mac[j];
					}
					open_bit_shares_for_plaintext_input_send[j].push_back(mbit);
				}
			}
		}

		vector<vector<BitWithMac>> open_bit_shares_for_plaintext_input_recv;
		open_bit_shares_for_plaintext_input_recv.resize(nP + 1);
		for(int j = 1; j <= nP; j++) {
			open_bit_shares_for_plaintext_input_recv[j].resize(len);
		}

		/*
		 * exchange the opening of the input mask
		 */
		vector<future<void>> res;
		for (int i = 1; i <= nP; ++i) {
			for (int j = 1; j <= nP; ++j) {
				if ((i < j) and (i == party or j == party)) {
					int party2 = i + j - party;

					res.push_back(pool->enqueue([this, &open_bit_shares_for_plaintext_input_recv, party2]() {
						io->recv_data(party2, open_bit_shares_for_plaintext_input_recv[party2].data(), sizeof(BitWithMac) * len);
						io->flush(party2);
					}));
					res.push_back(pool->enqueue([this, &open_bit_shares_for_plaintext_input_send, party2]() {
						io->send_data(party2, open_bit_shares_for_plaintext_input_send[party2].data(), sizeof(BitWithMac) * len);
						io->flush(party2);
					}));
				}
			}
		}
		joinNclean(res);

		/*
		 * verify the input mask
		 */
		vector<future<bool>> res_check;
		for (int j = 1; j <= nP; ++j) {
			if(j != party) {
				res_check.push_back(pool->enqueue([this, &input_mask, &open_bit_shares_for_plaintext_input_recv, j]() {
					for (int i = 0; i < len; i++) {
						if (party_assignment[i] == party) {
							block supposed_mac = Delta & select_mask[open_bit_shares_for_plaintext_input_recv[j][i].bit_share? 1 : 0];
							supposed_mac ^= input_mask[i].key[j];

							block provided_mac = open_bit_shares_for_plaintext_input_recv[j][i].mac;

							if(!cmpBlock(&supposed_mac, &provided_mac, 1)) {
								return true;
							}
						}
					}
					return false;
				}));
			}
		}
		if(joinNcleanCheat(res_check)) error("cheat in FlexIn's plaintext input mask!");

		/*
		 * broadcast the masked input
		 */
		vector<char> masked_input_sent;	// use char instead of bool because bools seem to fail for "data()"
		vector<vector<char>> masked_input_recv;
		masked_input_sent.resize(len);
		masked_input_recv.resize(nP + 1);
		for(int j = 1; j <= nP; j++) {
			masked_input_recv[j].resize(len);
		}

		for(int i = 0; i < len; i++) {
			if(party_assignment[i] == party) {
				masked_input_sent[i] = plaintext_assignment[i] ^ input_mask[i].bit_share;
				for(int j = 1; j <= nP; j++) {
					if(j != party) {
						masked_input_sent[i] = masked_input_sent[i] ^ open_bit_shares_for_plaintext_input_recv[j][i].bit_share;
					}
				}
			}
		}

		for (int i = 1; i <= nP; ++i) {
			for (int j = 1; j <= nP; ++j) {
				if ((i < j) and (i == party or j == party)) {
					int party2 = i + j - party;

					res.push_back(pool->enqueue([this, &masked_input_recv, party2]() {
						io->recv_data(party2, masked_input_recv[party2].data(), sizeof(char) * len);
						io->flush(party2);
					}));
					res.push_back(pool->enqueue([this, &masked_input_sent, party2]() {
						io->send_data(party2, masked_input_sent.data(), sizeof(char) * len);
						io->flush(party2);
					}));
				}
			}
		}
		joinNclean(res);

		vector<bool> masked_input;
		masked_input.resize(len);
		for(int i = 0; i < len; i++) {
			if(party_assignment[i] > 0) {
				int this_party = party_assignment[i];
				if(this_party == party) {
					masked_input[i] = masked_input_sent[i];
				} else {
					masked_input[i] = masked_input_recv[this_party][i];
				}
			}
		}

		/*
		 * secondly, handle the case party_assignment[] == -1
		 */

		/*
		 * Compute the authenticated bit share to the new circuit
		 * by XOR-ing with the input mask
		 */
		vector<AuthBitShare<nP>> authenticated_bitshares_new_circuit;
		for(int i = 0; i < len; i++) {
			AuthBitShare<nP> new_entry;
			memset(&new_entry, 0, sizeof(AuthBitShare<nP>));
			if(party_assignment[i] == -1) {
				new_entry.bit_share = authenticated_share_assignment[i].bit_share ^ input_mask[i].bit_share;
				for (int j = 1; j <= nP; j++) {
					new_entry.key[j] = authenticated_share_assignment[i].key[j] ^ input_mask[i].key[j];
					new_entry.mac[j] = authenticated_share_assignment[i].mac[j] ^ input_mask[i].mac[j];
				}
			}
			authenticated_bitshares_new_circuit.emplace_back(new_entry);
		}
		
		//print_block(Delta);
		
		/*
		cout << "Debug the authenticated input" << endl;
		for(int i = 0; i < 10; i++){
			if(party_assignment[i] == -1) {
				cout << "index: " << i << ", value: " << authenticated_share_assignment[i].bit_share << endl;
				cout << "mac: " << endl;
				for(int j = 1; j <= nP; j++) {
					if(j != party) {
						//cout << j << ": ";
						print_block(authenticated_share_assignment[i].mac[j]);
					}
				}
				cout << "key: " << endl;
				for(int j = 1; j <= nP; j++) {
					if(j != party) {
						//cout << j << ": ";
						print_block(authenticated_share_assignment[i].key[j]);
						print_block(authenticated_share_assignment[i].key[j] ^ Delta);
					}
				}
				cout << "=============" << endl;
			}
		}
		 */

		/*
		 * Opening the authenticated shares
		 */
		vector<vector<BitWithMac>> open_bit_shares_for_authenticated_bits_send;
		open_bit_shares_for_authenticated_bits_send.resize(nP + 1);

		for(int j = 1; j <= nP; j++) {
			if(j != party) {
				for(int i = 0; i < len; i++) {
					BitWithMac mbit{};
					if(party_assignment[i] == -1) {
						mbit.bit_share = authenticated_bitshares_new_circuit[i].bit_share;
						mbit.mac = authenticated_bitshares_new_circuit[i].mac[j];
					}
					open_bit_shares_for_authenticated_bits_send[j].push_back(mbit);
				}
			}
		}

		vector<vector<BitWithMac>> open_bit_shares_for_authenticated_bits_recv;
		open_bit_shares_for_authenticated_bits_recv.resize(nP + 1);
		for(int j = 1; j <= nP; j++) {
			open_bit_shares_for_authenticated_bits_recv[j].resize(len);
		}

		for (int i = 1; i <= nP; ++i) {
			for (int j = 1; j <= nP; ++j) {
				if ((i < j) and (i == party or j == party)) {
					int party2 = i + j - party;

					res.push_back(pool->enqueue([this, &open_bit_shares_for_authenticated_bits_recv, party2]() {
						io->recv_data(party2, open_bit_shares_for_authenticated_bits_recv[party2].data(), sizeof(BitWithMac) * len);
						io->flush(party2);
					}));
					res.push_back(pool->enqueue([this, &open_bit_shares_for_authenticated_bits_send, party2]() {
						io->send_data(party2, open_bit_shares_for_authenticated_bits_send[party2].data(), sizeof(BitWithMac) * len);
						io->flush(party2);
					}));
				}
			}
		}
		joinNclean(res);

		/*
		 * verify the input mask shares
		 */
		for (int j = 1; j <= nP; ++j) {
			if(j != party) {
				res_check.push_back(pool->enqueue([this, &authenticated_bitshares_new_circuit, &open_bit_shares_for_authenticated_bits_recv, j]() {
					for (int i = 0; i < len; i++) {
						if (party_assignment[i] == -1) {
							block supposed_mac = Delta & select_mask[open_bit_shares_for_authenticated_bits_recv[j][i].bit_share? 1 : 0];
							supposed_mac ^= authenticated_bitshares_new_circuit[i].key[j];

							block provided_mac = open_bit_shares_for_authenticated_bits_recv[j][i].mac;

							if(!cmpBlock(&supposed_mac, &provided_mac, 1)) {
								return true;
							}
						}
					}
					return false;
				}));
			}
		}
		if(joinNcleanCheat(res_check)) error("cheat in FlexIn's authenticated share input mask!");

		/*
		 * Reconstruct the authenticated shares
		 */
		for(int i = 0; i < len; i++) {
			if(party_assignment[i] == -1) {
				masked_input[i] = authenticated_bitshares_new_circuit[i].bit_share;
				for(int j = 1; j <= nP; j++) {
					if(j != party) {
						masked_input[i] = masked_input[i] ^ open_bit_shares_for_authenticated_bits_recv[j][i].bit_share;
					}
				}
			}
		}

		/*
		 * thirdly, handle the case party_assignment[] = -2
		 */

		/*
		 * Collect the masked input shares for un-authenticated bits
		 */
		vector<char> open_bit_shares_for_unauthenticated_bits_send;
		open_bit_shares_for_unauthenticated_bits_send.resize(len);

		for(int i = 0; i < len; i++) {
			if(party_assignment[i] == -2) {
				open_bit_shares_for_unauthenticated_bits_send[i] = (plaintext_assignment[i] ^ input_mask[i].bit_share)? 1 : 0;
			}
		}

		vector<vector<char>> open_bit_shares_for_unauthenticated_bits_recv;
		open_bit_shares_for_unauthenticated_bits_recv.resize(nP + 1);
		for(int j = 1; j <= nP; j++) {
			open_bit_shares_for_unauthenticated_bits_recv[j].resize(len);
		}

		for (int i = 1; i <= nP; ++i) {
			for (int j = 1; j <= nP; ++j) {
				if ((i < j) and (i == party or j == party)) {
					int party2 = i + j - party;

					res.push_back(pool->enqueue([this, &open_bit_shares_for_unauthenticated_bits_recv, party2]() {
						io->recv_data(party2, open_bit_shares_for_unauthenticated_bits_recv[party2].data(), sizeof(char) * len);
						io->flush(party2);
					}));
					res.push_back(pool->enqueue([this, &open_bit_shares_for_unauthenticated_bits_send, party2]() {
						io->send_data(party2, open_bit_shares_for_unauthenticated_bits_send.data(), sizeof(char) * len);
						io->flush(party2);
					}));
				}
			}
		}
		joinNclean(res);

		/*
		 * update the array of masked_input accordingly
		 */
		for(int i = 0; i < len; i++) {
			if(party_assignment[i] == -2) {
				masked_input[i] = open_bit_shares_for_unauthenticated_bits_send[i];
				for(int j = 1; j <= nP; j++) {
					if(j != party) {
						masked_input[i] = masked_input[i] ^ (open_bit_shares_for_unauthenticated_bits_recv[j][i] == 1);
					}
				}
			}
		}

		/*
		 * lastly, handle the case party_assignment[] = 0
		 */

		/*
		 * broadcast the input mask and its MAC
		 */
		vector<vector<BitWithMac>> open_bit_shares_for_public_input_send;
		open_bit_shares_for_public_input_send.resize(nP + 1);

		for(int j = 1; j <= nP; j++) {
			if(j != party) {
				for(int i = 0; i < len; i++) {
					BitWithMac mbit{};
					if(party_assignment[i] == 0) {
						mbit.bit_share = input_mask[i].bit_share;
						mbit.mac = input_mask[i].mac[j];
					}
					open_bit_shares_for_public_input_send[j].push_back(mbit);
				}
			}
		}

		vector<vector<BitWithMac>> open_bit_shares_for_public_input_recv;
		open_bit_shares_for_public_input_recv.resize(nP + 1);
		for(int j = 1; j <= nP; j++) {
			open_bit_shares_for_public_input_recv[j].resize(len);
		}

		/*
		 * exchange the opening of the input mask
		 */
		for (int i = 1; i <= nP; ++i) {
			for (int j = 1; j <= nP; ++j) {
				if ((i < j) and (i == party or j == party)) {
					int party2 = i + j - party;

					res.push_back(pool->enqueue([this, &open_bit_shares_for_public_input_recv, party2]() {
						io->recv_data(party2, open_bit_shares_for_public_input_recv[party2].data(), sizeof(BitWithMac) * len);
						io->flush(party2);
					}));
					res.push_back(pool->enqueue([this, &open_bit_shares_for_public_input_send, party2]() {
						io->send_data(party2, open_bit_shares_for_public_input_send[party2].data(), sizeof(BitWithMac) * len);
						io->flush(party2);
					}));
				}
			}
		}
		joinNclean(res);

		/*
		 * verify the input mask
		 */
		for (int j = 1; j <= nP; ++j) {
			if(j != party) {
				res_check.push_back(pool->enqueue([this, &input_mask, &open_bit_shares_for_public_input_recv, j]() {
					for (int i = 0; i < len; i++) {
						if (party_assignment[i] == 0) {
							block supposed_mac = Delta & select_mask[open_bit_shares_for_public_input_recv[j][i].bit_share? 1 : 0];
							supposed_mac ^= input_mask[i].key[j];

							block provided_mac = open_bit_shares_for_public_input_recv[j][i].mac;

							if(!cmpBlock(&supposed_mac, &provided_mac, 1)) {
								return true;
							}
						}
					}
					return false;
				}));
			}
		}
		if(joinNcleanCheat(res_check)) error("cheat in FlexIn's public input mask!");

		/*
		 * update the masked input
		 */
		for(int i = 0; i < len; i++) {
			if(party_assignment[i] == 0) {
				masked_input[i] = plaintext_assignment[i] ^ input_mask[i].bit_share;
				for(int j = 1; j <= nP; j++) {
					if(j != party) {
						masked_input[i] = masked_input[i] ^ open_bit_shares_for_public_input_recv[j][i].bit_share;
					}
				}
			}
		}

		/*
		 cout << "masked_input" << endl;
		for(int i = 0; i < len; i++) {
			cout << masked_input[i] << " ";
		}
		cout << endl;
		 */

		for(int i = 0; i < len; i++) {
			masked_input_ret[i] = masked_input[i];
		}
	}

	int get_length() {
		return len;
	}
};

template<int nP>
class FlexOut
{
public:

	int len{};
	int party{};

	bool cmpc_associated = false;
	ThreadPool* pool;
	bool *value;
	block *key[nP + 1];
	block *mac[nP + 1];
	block *eval_labels[nP + 1];
	NetIOMP<nP> * io;
	block Delta;
	block *labels;

	vector<int> party_assignment;
	// -1 represents an authenticated share,
	//  0 represents public output

	vector<bool> plaintext_results; // if `party` provides the value for this bit, the plaintext value is here
	vector<AuthBitShare<nP>> authenticated_share_results; // if this bit is from authenticated shares, the authenticated share is stored here

	FlexOut(int len, int party) {
		this->len = len;
		this->party = party;
		this->pool = pool;

		AuthBitShare<nP> empty_abit;
		memset(&empty_abit, 0, sizeof(AuthBitShare<nP>));

		party_assignment.resize(len, 0);
		plaintext_results.resize(len, false);
		authenticated_share_results.resize(len, empty_abit);
	}

	~FlexOut() {
		party_assignment.clear();
		plaintext_results.clear();
		authenticated_share_results.clear();
	}

	void associate_cmpc(ThreadPool *associated_pool, bool *associated_value, block *associated_mac[nP + 1], block *associated_key[nP + 1], block *associated_eval_labels[nP + 1], block *associated_labels, NetIOMP<nP> *associated_io, block associated_Delta) {
		this->cmpc_associated = true;
		this->pool = associated_pool;
		this->value = associated_value;
		this->labels = associated_labels;
		for (int j = 1; j <= nP; j++) {
			this->mac[j] = associated_mac[j];
			this->key[j] = associated_key[j];
		}
		if (party == ALICE){
			for (int j = 2; j <= nP; j++) {
				this->eval_labels[j] = associated_eval_labels[j];
			}
		}
		this->io = associated_io;
		this->Delta = associated_Delta;
	}

	void assign_party(int pos, int which_party) {
		party_assignment[pos] = which_party;
	}

	bool get_plaintext_bit(int pos) {
		assert(party_assignment[pos] == party || party_assignment[pos] == 0);
		return plaintext_results[pos];
	}
	
	AuthBitShare<nP> get_authenticated_bitshare(int pos) {
		assert(party_assignment[pos] == -1);
		return authenticated_share_results[pos];
	}

	int get_length() {
		return len;
	}

	void output(bool *masked_input_ret, int output_shift) {
		assert(cmpc_associated);

		/*
		 * Party 1 sends the labels of all the output wires out.
		 */
		vector<block> output_wire_label_recv;
		output_wire_label_recv.resize(len);

		vector<future<void>> res;

		if(party == ALICE) {
			vector<vector<block>> output_wire_label_send;
			output_wire_label_send.resize(nP + 1);

			for (int j = 2; j <= nP; j++) {
				output_wire_label_send[j].resize(len);
				for(int i = 0; i < len; i++) {
					output_wire_label_send[j][i] = eval_labels[j][output_shift + i];
				}
			}

			for(int j = 2; j <= nP; j++) {
				res.push_back(pool->enqueue([this, &output_wire_label_send, j]() {
					io->send_data(j, output_wire_label_send[j].data(), sizeof(block) * len);
					io->flush(j);
				}));
			}
			joinNclean(res);
		}else {
			io->recv_data(ALICE, output_wire_label_recv.data(), sizeof(block) * len);
			io->flush(ALICE);
		}

		/*
		 * Each party extracts x ^ r of each output wire
		 */
		vector<bool> masked_output;
		masked_output.resize(len);

		if(party == ALICE) {
			for(int i = 0; i < len; i++) {
				masked_output[i] = masked_input_ret[output_shift + i];
			}
		} else {
			for(int i = 0; i < len; i++) {
				block cur_label = output_wire_label_recv[i];
				block zero_label = labels[i + output_shift];
				block one_label = zero_label ^ Delta;

				if(cmpBlock(&cur_label, &zero_label, 1)) {
					masked_output[i] = false;
				} else if(cmpBlock(&cur_label, &one_label, 1)) {
					masked_output[i] = true;
				} else {
					error("Output label mismatched.\n");
				}
			}
		}

		/*
		 * Decide the broadcasting of the shares of r, as well as their MAC
		 */
		vector<vector<BitWithMac>> output_mask_send;
		output_mask_send.resize(nP + 1);
		for(int j = 1; j <= nP; j++) {
			BitWithMac empty_mbit{};
			output_mask_send[j].resize(len, empty_mbit);
		}

		for(int i = 0; i < len; i++) {
			if(party_assignment[i] == -1) {
				// do nothing, just update the share (for the first party) later
			} else if(party_assignment[i] == 0) {
				// public output, all parties receive the mbit
				for(int j = 1; j <= nP; j++){
					output_mask_send[j][i].bit_share = value[output_shift + i];
					output_mask_send[j][i].mac = mac[j][output_shift + i];
				}
			} else {
				// only one party is supposed to receive the mbit
				int cur_party = party_assignment[i];
				output_mask_send[cur_party][i].bit_share = value[output_shift + i];
				output_mask_send[cur_party][i].mac = mac[cur_party][output_shift + i];
			}
		}

		/*
		 * Exchange the output mask
		 */
		vector<vector<BitWithMac>> output_mask_recv;
		output_mask_recv.resize(nP + 1);
		for(int j = 1; j <= nP; j++) {
			output_mask_recv[j].resize(len);
		}

		for (int i = 1; i <= nP; ++i) {
			for (int j = 1; j <= nP; ++j) {
				if ((i < j) and (i == party or j == party)) {
					int party2 = i + j - party;

					res.push_back(pool->enqueue([this, &output_mask_recv, party2]() {
						io->recv_data(party2, output_mask_recv[party2].data(), sizeof(BitWithMac) * len);
						io->flush(party2);
					}));
					res.push_back(pool->enqueue([this, &output_mask_send, party2]() {
						io->send_data(party2, output_mask_send[party2].data(), sizeof(BitWithMac) * len);
						io->flush(party2);
					}));
				}
			}
		}
		joinNclean(res);

		/*
		 * Verify the output mask
		 */
		vector<future<bool>> res_check;
		for (int j = 1; j <= nP; ++j) {
			if(j != party) {
				res_check.push_back(pool->enqueue([this, &output_mask_recv, j, output_shift]() {
					for (int i = 0; i < len; i++) {
						if (party_assignment[i] == party || party_assignment[i] == 0) {
							block supposed_mac = Delta & select_mask[output_mask_recv[j][i].bit_share? 1 : 0];
							supposed_mac ^= key[j][output_shift + i];

							block provided_mac = output_mask_recv[j][i].mac;

							if(!cmpBlock(&supposed_mac, &provided_mac, 1)) {
								return true;
							}
						}
					}
					return false;
				}));
			}
		}
		if(joinNcleanCheat(res_check)) error("cheat in FlexOut's output mask!");

		/*
		 * Handle the case party_assignment[] = -1
		 */
		if(party == ALICE) {
			for(int i = 0; i < len; i++) {
				if(party_assignment[i] == -1) {
					authenticated_share_results[i].bit_share = value[output_shift + i] ^ masked_output[i];
					for(int j = 1; j <= nP; j++) {
						if(j != party) {
							authenticated_share_results[i].mac[j] = mac[j][output_shift + i];
							authenticated_share_results[i].key[j] = key[j][output_shift + i];
						}
					}
				}
			}
		} else {
			for(int i = 0; i < len; i++) {
				if(party_assignment[i] == -1) {
					authenticated_share_results[i].bit_share = value[output_shift + i];
					for(int j = 1; j <= nP; j++) {
						if(j != party) {
							authenticated_share_results[i].mac[j] = mac[j][output_shift + i];
							if(j == ALICE) {
								authenticated_share_results[i].key[j] =
										key[j][output_shift + i] ^ (Delta & select_mask[masked_output[i] ? 1 : 0]);
								// change the MAC key for the first party
							} else {
								authenticated_share_results[i].key[j] = key[j][output_shift + i];
							}
						}
					}
				}
			}
		}
		
		//print_block(Delta);
		
		/*
		cout << "Debug the authenticated output" << endl;
		for(int i = 0; i < 10; i++){
			if(party_assignment[i] == -1) {
				cout << "index: " << i << ", value: " << authenticated_share_results[i].bit_share << endl;
				cout << "mac: " << endl;
				for(int j = 1; j <= nP; j++) {
					if(j != party) {
						//cout << j << ": ";
						print_block(authenticated_share_results[i].mac[j]);
					}
				}
				cout << "key: " << endl;
				for(int j = 1; j <= nP; j++) {
					if(j != party) {
						//cout << j << ": ";
						print_block(authenticated_share_results[i].key[j]);
						print_block(authenticated_share_results[i].key[j] ^ Delta);
					}
				}
				cout << "=============" << endl;
			}
		}
		 */

		/*
		 * Handle the case party_assignment[] = 0 or == party
		 */
		for(int i = 0; i < len; i++) {
			if(party_assignment[i] == 0 || party_assignment[i] == party) {
				plaintext_results[i] = value[output_shift + i] ^ masked_output[i];
				for(int j = 1; j <= nP; j++) {
					if(j != party) {
						plaintext_results[i] = plaintext_results[i] ^ output_mask_recv[j][i].bit_share;
					}
				}
			}
		}
	}
};

#endif //EMP_AGMPC_FLEXIBLE_INPUT_OUTPUT_H
