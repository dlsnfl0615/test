#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <chrono>

using namespace std;


// ���̻� �迭 ����
vector<int> buildSuffixArray(const string& s) {
	int n = s.length();
	vector<int> sa(n); // ���̻��� ���� �ε����� �����ϴ� �迭
	vector<int> rank(n);
	/* rank[i]: ���̻� s[i:]�� �� k���ڿ� �ش��ϴ� ���� ����
	 * rank[i + k]: �� ���� k���� (s[i+k : i+2k])�� ���� ����
	 * (rank[i], rank[i + k])�� ���̻� s[i:]�� �� 2k����"�� �� �κ����� �ɰ��� ���� �Ͱ� ���� */
	vector<int> tmp(n);

	for (int i = 0; i < n; i++) {
		sa[i] = i; // [0,1,2,3,4,5,6]
		rank[i] = s[i]; // [98, 97, 110, 97, 110, 97, 36] = [b, a, n, a, n, a, $]
	}

	for (int k = 1; k < n; k *= 2) {
		auto cmp = [&](int i, int j) {
			if (rank[i] != rank[j]) return rank[i] < rank[j];
			// �� ���� ���ڳ��� ��
			int ri = i + k < n ? rank[i + k] : -1;
			int rj = j + k < n ? rank[j + k] : -1;
			return ri < rj;
			};
		sort(sa.begin(), sa.end(), cmp);

		tmp[sa[0]] = 0;
		for (int i = 1; i < n; i++) {
			tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
			/* ���� ���̻�� ������ ���� rank, �ٸ��� + 1 ����
			 * ���� ���̻簡 (1, 2)�̰� ���� ���̻簡 (1, 2)��� ������ ������ ����
			 * tmp[sa[i]] : sa[i]��°���� �����ϴ� ���̻��� ����. ex) temp[sa[2]] = 3�� ���� 2��°���� �����ϴ� ���̻��� ������ 3����� �� */
		}
		rank = tmp;
		// k = 4�� �� ���� 4���ھ��� ���ؾ��ϴµ� �̹� k = 2�� �� �� �� ������ ��� ������ ���س����Ƿ�
	}
	return sa;
}

// BWT ����
string buildBWT(const string& s, const vector<int>& sa) {
	string bwt;
	int n = s.length();
	for (int i = 0; i < n; i++) {
		bwt += (sa[i] == 0) ? s[n - 1] : s[sa[i] - 1];
	}
	return bwt;
}

// FM-Index Ŭ����
class FMIndex {
private:
	string bwt;
	vector<int> sa;
	map<char, int> first;
	vector<map<char, int>> occ;

public:
	FMIndex(const string& ref) {
		sa = buildSuffixArray(ref + '$');
		bwt = buildBWT(ref + '$', sa);

		map<char, int> char_count;
		for (char c : bwt) char_count[c]++;
		int sum = 0;
		for (auto& p : char_count) {
			first[p.first] = sum;
			sum += p.second;
		}

		occ.resize(bwt.length());
		map<char, int> running_count;
		for (size_t i = 0; i < bwt.length(); i++) {
			running_count[bwt[i]]++;
			occ[i] = running_count;
		}
	}

	// ��Ȯ�� ��Ī�� ���� ���� ���
	bool getRange(char c, int& top, int& bottom) {
		if (first.find(c) == first.end()) return false;
		top = first[c] + (top > 0 ? occ[top - 1][c] : 0);
		bottom = first[c] + occ[bottom][c] - 1;
		return top <= bottom;
	}

	// �ٻ� ��Ī (�ִ� k���� mismatch ���)
	void approxSearch(const string& pattern, int k, int pos, int top, int bottom, vector<int>& results) {
		if (pos < 0) {
			for (int i = top; i <= bottom; i++) {
				if (sa[i] < pattern.length()) continue; // '$' ����
				results.push_back(sa[i]);
			}
			return;
		}

		// ��Ȯ�� ��Ī �õ�
		int new_top = top, new_bottom = bottom;
		if (getRange(pattern[pos], new_top, new_bottom)) {
			approxSearch(pattern, k, pos - 1, new_top, new_bottom, results);
		}

		// k�� ���������� mismatch �õ�
		if (k > 0) {
			for (char c : {'A', 'C', 'G', 'T'}) { // DNA ���� ����
				if (c == pattern[pos]) continue;
				new_top = top;
				new_bottom = bottom;
				if (getRange(c, new_top, new_bottom)) {
					approxSearch(pattern, k - 1, pos - 1, new_top, new_bottom, results);
				}
			}
		}
	}

	vector<int> searchWithMismatch(const string& pattern, int k) {
		vector<int> results;
		approxSearch(pattern, k, pattern.length() - 1, 0, bwt.length() - 1, results);
		sort(results.begin(), results.end());
		return results;
	}
};

// ���� �б� �Լ�
string readReference(const string& filename) {
	ifstream file(filename);
	string ref;
	getline(file, ref);
	file.close();
	return ref;
}

vector<string> readPatterns(const string& filename) {
	ifstream file(filename);
	vector<string> patterns;
	string line;
	while (getline(file, line)) {
		patterns.push_back(line);
	}
	file.close();
	return patterns;
}

vector<int> readGroundTruth(const string& filename) {
	ifstream file(filename);
	vector<int> truth;
	int pos;
	while (file >> pos) {
		truth.push_back(pos);
	}
	file.close();
	return truth;
}

int main() {
	auto start = chrono::high_resolution_clock::now();

	// ���� �б�
	string ref = readReference("reference_1M.txt");
	vector<string> patterns = readPatterns("mammoth_reads_10K.txt");
	vector<int> ground_truth = readGroundTruth("ground_truth_10K.txt");

	// FM-Index ����
	FMIndex fm(ref);

	// �ִ� mismatch ��� ��
	int k = 2; // �ʿ�� ���� ����

	// ��� ��� ����
	vector<string> outputs;
	ofstream out("approx_search_results.txt");

	// ���� �˻� �� ����
	int correct = 0;
	for (size_t i = 0; i < patterns.size(); i++) {
		vector<int> positions = fm.searchWithMismatch(patterns[i], k);
		bool found = false;
		for (int pos : positions) {
			if (pos == ground_truth[i]) {
				found = true;
				break;
			}
		}
		if (found) correct++;
		outputs.push_back(("Pattern " + to_string(i + 1) + ": " + (found ? "Match at " + to_string(ground_truth[i]) : "No match")));
	}
	auto end = chrono::high_resolution_clock::now();

	cout << "Accuracy: " << (double)correct / patterns.size() * 100 << "%" << endl;
	for (const string& output : outputs) {
		out << output << endl;
	}
	out << "Accuracy: " << (double)correct / patterns.size() * 100 << "%" << endl;
	out.close();

	
	auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
	cout << "Execution time: " << duration.count() << " ms" << endl;

	return 0;
}